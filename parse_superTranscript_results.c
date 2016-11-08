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

//minimum number of bases of gap reported for query to be considered interesting
//(e.g. bases in assembly not in genome). Indel size
static const int qGapMin=7;
//minimum number of bases at the edge of the contig not found in the genome
//best to make this > smallest size that can be aligned with blat
static const int qEdgeMin=30;

// the real stuff starts here.
int main(int argc, char **argv){

  if(argc!=4){
    cout << "Wrong number of arguments" << endl;
    exit(1);
  }

  ifstream file;
  //Open the exons file
  file.open(argv[1]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  // reading the known junction list file
  string line;
  string column;
  string chrom;
  vector<g_interval> known_junc;
  for(int skip=1; skip!=0 && getline(file,line) ; skip--); //skip the first line
  while(getline(file,line) ){
    istringstream line_stream(line);
    line_stream >> column; //transcript ID (not used)
    line_stream >> chrom;
    line_stream >> column;
    vector<int> junc_start = get_vector_from_list(column);
    line_stream >> column;
    vector<int> junc_end = get_vector_from_list(column);
    for(int i=0; i < (junc_end.size()-1) ; i++){
      g_interval new_junc;
      new_junc.chrom=chrom;
      new_junc.start=junc_end.at(i);
      new_junc.end=junc_start.at(i+1);
      known_junc.push_back(new_junc);
    }
  }
  file.close();
  
  /***************************
   ** now read the blat file
   ***************************/
  file.open(argv[2]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[2] << endl;
    exit(1);
  }
  vector<string> non_interesting_ids;
  vector<string> all_ids;
  for(int skip=5; skip!=0 && getline(file,line) ; skip--); //skip the first line
  while(getline(file,line) ){
    istringstream line_stream(line);
    vector<string> buffer;
    string value;
    while( line_stream >> value){ 
      buffer.push_back(value);
    };
    all_ids.push_back(buffer[9]);
    bool novel_seq = atoi(buffer[5].c_str()) >= qGapMin ;
    int contig_length=atoi(buffer[10].c_str()) ;
    int contig_start=atoi(buffer[11].c_str()) ;
    int contig_end=atoi(buffer[12].c_str()) ;
    bool edge_missing =
      ((contig_length - contig_end ) >= qEdgeMin ) ||
      (contig_start >= qEdgeMin ) ;
    // code where we work out the junctions...
    string chrom=buffer[13];
    vector<int> junc_start=get_vector_from_list(buffer[20]);
    vector<int> junc_size=get_vector_from_list(buffer[18]);
    string pos;
    //now check if the junctions are known
    bool is_junc_known=true;
    for(int i=0; i < (junc_size.size()-1) ; i++){
      g_interval this_junc;
      this_junc.chrom=chrom;
      this_junc.start=junc_start.at(i)+junc_size.at(i);
      this_junc.end=junc_start.at(i+1);
      if((this_junc.end - this_junc.start) >= qGapMin){
	is_junc_known=find(known_junc.begin(),known_junc.end(),this_junc) != known_junc.end();
	if(!is_junc_known){
	  //	  cout << buffer[9] << " " << chrom << " " << this_junc.start << " " << this_junc.end << endl;
	  break;
	}
      }
    }
    if(!novel_seq & !edge_missing & is_junc_known) //add id to the exclude list
      non_interesting_ids.push_back(buffer[9]);
  }
  file.close(); 
  sort_vector(all_ids) ;
  sort_vector(non_interesting_ids) ;
  
  vector<string> interesting_ids;
  set_difference(all_ids.begin(), all_ids.end(),
		 non_interesting_ids.begin(), non_interesting_ids.end(),
		 back_inserter( interesting_ids )
		 );

  file.open(argv[3]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[3] << endl;
    exit(1);
  }
  bool report=false;
  while ( getline (file,line) ){
    int start=line.find(">")+1;
    if(start==1){ //if this is the ID line...
      int end=line.find_first_of("\t\n ")-1;
      string id=line.substr(start,end);
      bool in_list=find(interesting_ids.begin(),interesting_ids.end(),id) != interesting_ids.end() ;
      in_list ? report=true : report=false ;
    }
    if(report) cout << line << endl; //output
    }
  return(0);
}

