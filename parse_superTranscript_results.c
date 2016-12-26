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
 * A list of novel splice junctions
 * A list of novel blocks
 * A list of the input contigs
 *
 * The program will output:
 * a list of filtered events
 **/

struct event {
  int start, end;
  int depth;
  string type;
} ;

static const int read_length=100;
static const double junc_overhand_frac=0.76; // read_length - 12 bases * 2 ends / read_length
static const int min_abs_depth=5;
static const double min_rel_depth=0.4;

// the real stuff starts here.
int main(int argc, char **argv){

  if(argc!=4){
    cout << "Wrong number of arguments" << endl;
    cout << "Usage: parse_superTranscript_results <novel.junctions> <novel.blocks> <all.groupings>" << endl;
    exit(1);
  }

  ifstream file;
  string line;
  string gene;
  string column;
  map< string , vector<event> > novel_events ;
  map< string, int > max_depth ;

  //Open the exons file
  file.open(argv[1]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  // reading the novel junction list file

  while(getline(file,line) ){
    istringstream line_stream(line);
    event new_event;
    line_stream >> gene;
    line_stream >> new_event.start;
    line_stream >> new_event.end;
    for(int skip=3; skip!=0 && (line_stream >> column) ; skip--); //skip some columns
    line_stream >> new_event.depth;
    new_event.type="deletion" ;
    novel_events[gene].push_back(new_event);
  }
  file.close();

  // reading the novel blocks list file
  file.open(argv[2]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[2] << endl;
    exit(1);
  }
  string block_type;
  while(getline(file,line) ){
    istringstream line_stream(line);
    event new_event;
    line_stream >> gene;
    line_stream >> new_event.start;
    line_stream >> new_event.end;
    line_stream >> block_type;
    line_stream >> new_event.depth;
    //fix the depth
    new_event.depth = new_event.depth / (new_event.end - new_event.start ) ;
    new_event.type="insertion" ;
    if(block_type=="Assembly")
      novel_events[gene].push_back(new_event);
    if(max_depth[gene]<new_event.depth) max_depth[gene]=new_event.depth;
  }
  file.close();

  //read the cluster file, so we can report which contigs went into the lace assembly
  file.open(argv[3]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[3] << endl;
    exit(1);
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

  //loop through the events and print out
  map<string,vector<event> >::iterator gItr=novel_events.begin();
  map<string,vector<event> >::iterator gItrEnd=novel_events.end();
  for(; gItr!=gItrEnd ; gItr++){
    vector<event> these_events=gItr->second;
    for(int e=0; e < these_events.size(); e++){
      bool has_abs_depth = these_events.at(e).depth >= min_abs_depth ;
      double rel_depth = these_events.at(e).depth / (double) max_depth[gItr->first] ;
      if(these_events.at(e).type=="deletion")
	rel_depth=(rel_depth/junc_overhand_frac);
      bool has_rel_depth = rel_depth >= min_rel_depth ;
      if(has_abs_depth && has_rel_depth){
	cout << gItr->first << "\t" 
	     << these_events.at(e).type << "\t"
	     << these_events.at(e).end - these_events.at(e).start
	     << these_events.at(e).depth << "\t"
	     << rel_depth << "\t"
	     << these_events.at(e).start << "\t" 
	     << these_events.at(e).end << "\t" ;
	//loop over the contigs assembled for this gene
	vector<string> conts = contigs[ gItr->first ];
	for( int c=0 ; c < conts.size() ; c++)
	  cout << conts.at(c) << "," ;
	cout << endl;
      }
    }
  }
  
  return(0);
}

