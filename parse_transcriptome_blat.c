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
  string line, temp, trans_id, gene_id;
  map<string,string> id_map;
  //  cout << "Reading reference file: " << argv[1] << endl;
  for(int skip=1; skip!=0 && getline(file,line) ; skip--); //skip the first line
  while(getline(file,line) ){
    istringstream line_stream(line);
    line_stream >> trans_id; //transcript ID 
    for(int i=0; i<3 ; i++) line_stream >> temp;
    line_stream >> gene_id;
    id_map[trans_id]=gene_id;
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
  map<string,vector<string> > gene_match;
  //cout << "Reading blat file: " << argv[2] << endl;
  for(int skip=5; skip!=0 && getline(file,line) ; skip--); //skip the first line
  //int i=0;
  while(getline(file,line) ){
    // i++;
    //    if( (i % 100) == 0) cout << "Alignment "<< i << endl;
    istringstream line_stream(line);
    vector<string> buffer;
    string value;
    while( line_stream >> value){ 
      buffer.push_back(value);
    };
    all_ids.push_back(buffer[9]);
    bool qGap = atoi(buffer[5].c_str()) >= qGapMin ;
    bool tGap = atoi(buffer[7].c_str()) >= qGapMin ;
    int contig_length=atoi(buffer[10].c_str()) ;
    int contig_start=atoi(buffer[11].c_str()) ;
    int contig_end=atoi(buffer[12].c_str()) ;
    bool edge_missing =
      ((contig_length - contig_end ) >= qEdgeMin ) ||
      (contig_start >= qEdgeMin ) ;
    if(!qGap & !tGap & !edge_missing){ //add id to the exclude list
      non_interesting_ids.push_back(buffer[9]);
    }
    string gene=id_map[buffer[13]];
    if(gene!="") gene_match[buffer[9]].push_back(gene);
    //    cout << buffer[9] << " " << buffer[13] << " " << id_map[buffer[13]] << endl;
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
  map<string,vector<string> > gene_cluster ;
  // cout << "Reporting interesting contigs" << endl;
  while ( getline (file,line) ){
    int start=line.find(">")+1;
    if(start==1){ //if this is the ID line...
      int end=line.find_first_of("\t\n ")-1;
      string id=line.substr(start,end);
      bool in_list=find(interesting_ids.begin(),interesting_ids.end(),id) != interesting_ids.end() ;
      if(in_list){ 
	report=true ;
	sort_vector(gene_match[id]) ;
	if(gene_match[id].size()>=1){
	  stringstream cluster;
	  cluster << gene_match[id].at(0);
	  for(int i=1; i<gene_match[id].size() ; i++) cluster << "|" << gene_match[id].at(i) ;
	  cout << id << "\t" << cluster.str() << endl;
	  for(int i=0; i<gene_match[id].size() ; i++) gene_cluster[gene_match[id].at(i)].push_back(cluster.str());
	}
      } else {
	report=false;
      }
    }
    //    if(report) cout << line << endl; //output
    }
  map<string,vector<string> >::iterator clustItr=gene_cluster.begin();
  map<string,vector<string> >::iterator clustItrEnd=gene_cluster.end();
  for( ; clustItr!=clustItrEnd ; clustItr++){
    vector<string> clusters = clustItr->second;
    sort_vector(clusters);
    for(int i=0; i<clusters.size(); i++){
      cout << clustItr->first << "\t" << clusters.at(i) << endl;
    } 
  }
  return(0);
}

