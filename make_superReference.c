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
	 !(y1<x2) && !(y2<x1) &&
	 interval.at(i).strand==interval.at(j).strand){
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

  if(argc!=4){
    cout << "Wrong number of arguments" << endl;
    cout << "Usage: " << endl;
    cout << "   make_superReference <in_uscs_table_file> <out_genome_gtf_flattened> <out_super_transcript_gtf>" << endl;
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
  cout << "Reading the input file..." << endl;
  string line;
  string column;
  string chrom;
  string symbol;
  string strand;
  map<string,vector<g_interval> > exons;
  for(int skip=1; skip!=0 && getline(file,line) ; skip--); //skip the first line
  while(getline(file,line) ){
    istringstream line_stream(line);
    line_stream >> column; //transcript ID (not used)
    line_stream >> chrom;
    line_stream >> strand;
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
      new_exon.strand=strand;
      exons[symbol].push_back(new_exon);
    }
  }
  file.close();

  const static int MAX_INTRON = 1000000;

  //loop over exons
  map<string,vector<g_interval> >::iterator exonItr=exons.begin();
  map<string,vector<g_interval> >::iterator exonItrEnd=exons.end();
  //Open the output file
  cout << "Merging exon intervals to make a flat gff file ..." << endl;
  ofstream flat_file, ann_file;
  static int chain_id=1;
  flat_file.open(argv[2]);
  ann_file.open(argv[3]);
  for( ; exonItr!=exonItrEnd ; exonItr++){ //looping over the genes
    vector<g_interval> exonInterval = merge_intervals(exonItr->second); //merging exons
    sort(exonInterval.begin(),exonInterval.end(),g_interval_compare);
    int last_start, last_end=0;
    vector<int> genome_starts,genome_ends,block_sizes;
    string chrom=exonInterval.at(0).chrom;
    string strand=exonInterval.at(0).strand;
    int st_length=0;
    for(int i=0; i<exonInterval.size(); i++){ //loop over the exons in the gene
      string this_chrom=exonInterval.at(i).chrom;
      int this_start=exonInterval.at(i).start;
      int this_end=exonInterval.at(i).end;
      string this_strand=exonInterval.at(i).strand;
      if(last_end==0 
	 || ( chrom==this_chrom && abs(this_start-last_end)<MAX_INTRON && strand==this_strand)){
	genome_starts.push_back(this_start);
	genome_ends.push_back(this_end);
	st_length+=this_end-this_start+1;
	block_sizes.push_back(this_end-this_start+1);
	last_start=this_start; last_end=this_end; 
      }
    }
    //now loop again over the filtered exons
    ann_file << "chain\t1\t" << chrom << "\t10000000000\t" 
	     << "\t+\t" << genome_starts.at(0)-1 
	     << "\t" << genome_ends.at(genome_ends.size()-1) 
	     << "\t" << exonItr->first << "\t" << st_length 
	     << "\t" << strand << "\t0\t" << st_length 
	     << "\t" << chain_id << endl;    
    for(int i=0; i<genome_starts.size(); i++){ 
      flat_file << chrom << "\t" ;
      flat_file << "superT\texon\t";
      flat_file << genome_starts.at(i) << "\t" ;
      flat_file << genome_ends.at(i) << "\t" ;
      flat_file << ".\t"<< strand << "\t.\t";
      flat_file << exonItr->first << endl;
      int i_max=genome_starts.size()-1;
      ann_file << genome_ends.at(i) - genome_starts.at(i) + 1 ;
      if(i<i_max)
	ann_file << "\t" << genome_starts.at(i+1) - genome_ends.at(i) - 1 << "\t0";
      ann_file << endl;
    }
    chain_id++;
  }
  flat_file.close();
  ann_file.close();

  cout << "All done" << endl;
}
