/**
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 26 October 2016
 **
 ** e.g. gffread temp.gff -g /group/bioi1/shared/genomes/hg38/fasta/hg38.fa -w temp.fasta
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

vector<g_interval> merge_intervals( vector<g_interval> interval){
  //loop over each interval
  vector<g_interval> result;
  if(interval.size()<2) return interval;
  for(int i=0; i<(interval.size()-1); i++){
    bool merged=false;
    for(int j=(i+1); (j<interval.size()) && !merged; j++){
      //check if the two intervals overlap
      int x1=interval.at(i).start;
      int y1=interval.at(i).end;
      int x2=interval.at(j).start;
      int y2=interval.at(j).end;
      if((interval.at(i).chrom==interval.at(j).chrom) &&
	 !(y1<x2) && !(y2<x1)){
	interval.at(j).start=min(x1,x2);
	interval.at(j).end=max(y1,y2);
	merged=true;
      }
    }
    if(!merged) result.push_back(interval.at(i));
  }
  result.push_back(interval.at(interval.size()-1));
  return result;
}

// the real stuff starts here.
int main(int argc, char **argv){

  if(argc!=2){
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
  string symbol;
  map<string,vector<g_interval> > exons;
  for(int skip=1; skip!=0 && getline(file,line) ; skip--); //skip the first line
  while(getline(file,line) ){
    istringstream line_stream(line);
    line_stream >> column; //transcript ID (not used)
    line_stream >> chrom;
    line_stream >> column;
    vector<int> starts = get_vector_from_list(column);
    line_stream >> column;
    vector<int> ends = get_vector_from_list(column);
    line_stream >> symbol;
    for(int i=0; i < starts.size() ; i++){
      g_interval new_exon;
      new_exon.chrom=chrom;
      new_exon.start=starts.at(i)+1;
      new_exon.end=ends.at(i);
      exons[symbol].push_back(new_exon);
    }
  }
  const static int MAX_INTRON = 1000000;

  //loop over exons
  map<string,vector<g_interval> >::iterator exonItr=exons.begin();
  map<string,vector<g_interval> >::iterator exonItrEnd=exons.end();
  for( ; exonItr!=exonItrEnd ; exonItr++){
    vector<g_interval> exonInterval = merge_intervals(exonItr->second);
    string last_chrom;
    int last_start, last_end=0;
    for(int i=0; i<exonInterval.size(); i++){
      string this_chrom=exonInterval.at(i).chrom;
      int this_start=exonInterval.at(i).start;
      int this_end=exonInterval.at(i).end;
      if(last_end==0 || ( last_chrom==this_chrom && abs(this_start-last_end)<MAX_INTRON)){
	cout << exonInterval.at(i).chrom << "\t" ;
	cout << "superT\texon\t";
	cout << exonInterval.at(i).start << "\t" ;
	cout << exonInterval.at(i).end << "\t" ;
	cout << ".\t+\t.\t";
	cout << exonItr->first << endl;
	last_start=this_start; last_end=this_end; last_chrom=this_chrom;
      }
    }
  }
  file.close();
  
}
