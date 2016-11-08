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
 * Take a psl file from blat and report the gtf * 
 **/
 

// the real stuff starts here.
int main(int argc, char **argv){
  if(argc!=2){
    cout << "Wrong number of arguments" << endl;
    cout << endl;
    cout << "Usage: psl3gtf <psl>" << endl;
    exit(1);
  }

  ifstream file;
  //Open the exons file
  file.open(argv[1]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  // reading the blat table
  string line;
  string column;
  string gene_id;
  string trans_id;
  string strand;
  for(int line_skip=5; line_skip!=0 && getline(file,line) ; line_skip--); //skip the first line
  while(getline(file,line) ){
    istringstream line_stream(line);
    for(int col_skip=8; col_skip!=0 && (line_stream >> column); col_skip--);
    line_stream >> strand;
    line_stream >> trans_id;
    for(int col_skip=3; col_skip!=0 && (line_stream >> column); col_skip--);
    line_stream >> gene_id;
    for(int col_skip=4; col_skip!=0 && (line_stream >> column); col_skip--);
    line_stream >> column;
    vector<int> junc_size = get_vector_from_list(column);
    line_stream >> column;
    line_stream >> column;
    vector<int> junc_start = get_vector_from_list(column);
    for(int i=0; i < junc_size.size()-1 ; i++){
      cout << gene_id << "\t"
	   << junc_start.at(i)+junc_size.at(i)+1 << "\t"
	   << junc_start.at(i+1) << "\t"
	   << endl;
    }
  }
  file.close();
}

