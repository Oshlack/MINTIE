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

using namespace std;

struct g_interval {
  string chrom;
  int start, end;
  int support;
  bool operator==(const g_interval& other)
  {
    return ((chrom==other.chrom) && (start==other.start) && (end==other.end));
  };
} ;

void sort_vector( vector<string> & a){ 
  sort(a.begin(), a.end());
  a.resize(distance(a.begin(),unique(a.begin(),a.end())));
};

vector<int> get_vector_from_list(string comm_list){
  vector<int> result;
  istringstream string_stream(comm_list);
  string pos;
  while(getline(string_stream,pos,',')){ 
    result.push_back(atoi(pos.c_str()));
  };
  return result;
};
