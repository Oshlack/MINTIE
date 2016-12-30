/**
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 26 October 2016
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "misc.h"

using namespace std;

/**
 * This program takes two files:
 * A list of novel splice junctions files
 * A list of novel blocks files
 * A list of the input contigs
 *
 * The program will output:
 * a list of filtered events
 **/

struct event {
  int start, end;
  int abs_depth;
  double rel_depth;
  string type;
  string sample;
} ;

// threasholds for filtering
static const double junc_overhand_frac=0.76; // read_length - 12 bases * 2 ends / read_length
static const int min_abs_depth=5;
static const double min_rel_depth=0.4;
static const int max_abs_depth_controls=5;
static const double max_rel_depth_controls=0.05;


map< string , vector<event> > get_files(const string block_file, const string junction_file){

  ifstream file;
  string line;
  string gene;
  string column;
  map< string , vector<event> > novel_events ;

  //Read the novel junction file
  file.open(junction_file.c_str());
  if(!(file.good())){
    cout << "Unable to open file " << junction_file << endl;
    exit(1);
  } else {
    //  cout << "Opening "<< junction_file << endl;
  }
  while(getline(file,line) ){
    istringstream line_stream(line);
    event new_event;
    line_stream >> gene;
    line_stream >> new_event.start;
    line_stream >> new_event.end;
    for(int skip=3; skip!=0 && (line_stream >> column) ; skip--); //skip some columns
    line_stream >> new_event.abs_depth;
    new_event.type="deletion" ;
    novel_events[gene].push_back(new_event);
  }
  file.close();

  //Read the novel block list file
  map< string, int > max_depth ;
  file.open(block_file.c_str());
  if(!(file.good())){
    cout << "Unable to open file " << block_file << endl;
    exit(1);
  } else {
    //   cout << "Opening "<< block_file << endl;
  }
  string block_type;
  while(getline(file,line) ){
    istringstream line_stream(line);
    event new_event;
    line_stream >> gene;
    line_stream >> new_event.start;
    line_stream >> new_event.end;
    line_stream >> block_type;
    line_stream >> new_event.abs_depth;
    new_event.type="insertion" ;
    if(block_type=="Assembly")
      novel_events[gene].push_back(new_event);
    if(max_depth[gene]<new_event.abs_depth) max_depth[gene]=new_event.abs_depth;
  }
  file.close();

  //Loop through the events and calculate the relative coverage 
  //(relative to coverage of the block with the maximum coverage)
  map<string,vector<event> >::iterator gItr=novel_events.begin();
  map<string,vector<event> >::iterator gItrEnd=novel_events.end();
  for(; gItr!=gItrEnd ; gItr++){
    vector<event> & these_events=gItr->second;
    for(int e=0; e < these_events.size(); e++){
      these_events.at(e).rel_depth = these_events.at(e).abs_depth / (double) max_depth[gItr->first] ;
      if(these_events.at(e).type=="deletion")
	these_events.at(e).rel_depth=(these_events.at(e).rel_depth/junc_overhand_frac);
    }
  }
  return novel_events;
}

// the real stuff starts here.
int main(int argc, char **argv){

  ifstream file;
  string line;
  string gene;
  string column;

  /**  if(argc<5){
    cout << "Wrong number of arguments" << endl;
    cout << "Usage: parse_superTranscript_results <novel.blocks> <novel.junctions> <all.groupings> <control files>" << endl;
    exit(1);
    }**/

  map< string , vector<event> > novel_events = get_files(argv[1],argv[2]);
  map< string , vector<event> > control_events;
  for(int a=4; a < argc ; a+=2){
    map< string , vector<event> > temp_events;
    temp_events = get_files(argv[a],argv[a+1]);
    //merge
    map<string,vector<event> >::iterator tItr=temp_events.begin();
    map<string,vector<event> >::iterator tItrEnd=temp_events.end();
    //Loop through the events
    for(; tItr!=tItrEnd ; tItr++){
      control_events[tItr->first].insert(control_events[tItr->first].end(), 
					 tItr->second.begin(),tItr->second.end() );
    }
  }

  //Filter the events..
  map< string , vector<event> > filtered_events;
  map<string,vector<event> >::iterator gItr=novel_events.begin();
  map<string,vector<event> >::iterator gItrEnd=novel_events.end();
  //Loop through the events and print out
  for(; gItr!=gItrEnd ; gItr++){
    vector<event> these_events=gItr->second;
    for(int e=0; e < these_events.size(); e++){
      // check that the depth is good
      if(these_events.at(e).abs_depth < min_abs_depth) break ;
      if(these_events.at(e).rel_depth < min_rel_depth) break ;
      // compare to the control samples
      int control_max_abs_depth=0;
      double control_max_rel_depth=0;
      vector<event> these_control_events=control_events[gItr->first];
      for(int c=0; c < these_control_events.size(); c++){
	if(these_events.at(e).end == these_control_events.at(c).end &&
	   these_events.at(e).start == these_control_events.at(c).start &&
	   these_events.at(e).type == these_control_events.at(c).type){
	  if( these_control_events.at(c).abs_depth> control_max_abs_depth) 
	    control_max_abs_depth=these_control_events.at(c).abs_depth;
	  if( these_control_events.at(c).rel_depth> control_max_rel_depth) 
	    control_max_rel_depth=these_control_events.at(c).rel_depth;
	}
      }
      if(control_max_abs_depth > max_abs_depth_controls){ break; }
      if(control_max_rel_depth > max_rel_depth_controls) break;
	filtered_events[gItr->first].push_back(these_events.at(e));
    }
  }
  
  //Read the cluster file, so we can report which contigs went into the lace assembly
  file.open(argv[3]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[3] << endl;
    exit(1);
  } else {
    //    cout << "Opening " << argv[3] << endl;
  }
  map < string, vector <string > > contigs;
  while(getline(file,line) ){
    istringstream line_stream(line);
    string cont;
    string gene;
    line_stream >> cont;
    line_stream >> gene;
    contigs[gene].push_back(cont);
  }
  file.close();

  gItr=filtered_events.begin();
  gItrEnd=filtered_events.end();
  //Loop through the events and print out
  for(; gItr!=gItrEnd ; gItr++){
    vector<event> these_events=gItr->second;
    for(int e=0; e < these_events.size(); e++){
	cout << gItr->first << "\t" 
	     << these_events.at(e).type << "\t"
	     << these_events.at(e).end - these_events.at(e).start << "\t"
	     << these_events.at(e).abs_depth << "\t"
	     << these_events.at(e).rel_depth << "\t"
	     << these_events.at(e).start << "\t" 
	     << these_events.at(e).end << "\t" ;
	//loop over the contigs assembled for this gene
	vector<string> conts = contigs[ gItr->first ];
	for( int c=0 ; c < conts.size() ; c++)
	  cout << conts.at(c) << "," ;
	cout << endl;
    }
  }
  
  return(0);
}

