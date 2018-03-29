#########################################################
# Author: Marek Cmero
# Create a bed annotation from a supertranscript's blat
# results (against the genome) in order to visualise genes
# and exons
#########################################################

import argparse
import pandas as pd
import numpy as np
import os
from intervaltree import Interval, IntervalTree
from random import randint

parser = argparse.ArgumentParser()
#parser.add_argument(dest='blat',
#                    help='''blat hits of the supertranscript against the genome''')
parser.add_argument(dest='gtf',
                    help='''gene reference in GTF format''')

args        = parser.parse_args()
#blat_file   = args.blat
gtf_file    = args.gtf
#out_prefix  = os.path.basename(blat_file)

gtf = pd.read_csv(gtf_file, header=None, sep='\t', comment='#')
#blat = pd.read_csv(blat_file, header=None, sep='\t', skiprows=5)
#blat_header = pd.read_csv(blat_file, header=None, sep='\t', nrows=2, skiprows=2)
#blat_header = blat_header.applymap(lambda x: '' if str(x)=='nan' else str(x).strip())
#blat_header = blat_header.apply(lambda x: ''.join(x))
#blat.columns = blat_header

#genes = gtf[gtf[2]=='gene']
#exons = gtf[gtf[2]=='exon']

st_ranges = {}
chroms, starts, ends = gtf[0].values, gtf[3].values-1, gtf[4].values
exons = [int(info.split(';')[2].strip().split(' ')[1]) for info in gtf[8]]
for i in range(len(chroms)):
   
    RGB = ','.join([str(randint(0,255)) for i in range(3)])
    line = '\t'.join([chroms[i], str(starts[i]), str(ends[i]), 'Block%s' % exons[i], '0', '.', str(starts[i]), str(ends[i]), RGB])
    print(line)    

#gene_ranges = {}
#chroms, starts, ends = genes[0].values, genes[3].values, genes[4].values
#gene_names = [info.split(';')[3].strip().split(' ')[1].replace('"','') for info in genes[8]]
#for i in range(len(chroms)):
#    if chroms[i] not in gene_ranges:
#        gene_ranges[chroms[i]] = IntervalTree()
#    gene_ranges[chroms[i]][starts[i]:ends[i]] = gene_names[i] 
#
#exon_ranges = {}
#chroms, starts, ends = exons[0].values, exons[3].values, exons[4].values
##exon_names = [(info.split(';')[9].strip().split(' ')[1].replace('"',''), info.split(';')[4].strip().split(' ')[1].replace('"',''))  for info in exons[8]]
#exon_names = [info.split(';')[9].strip().split(' ')[1].replace('"','') for info in exons[8]]
#for i in range(len(exon_names)):
#    if starts[i] >= ends[i]:
#        continue
#    if chroms[i] not in exon_ranges:
#        exon_ranges[chroms[i]] = IntervalTree()
#    exon_id = exon_names[i]
#    exon_ranges[chroms[i]][starts[i]:ends[i]] = exon_id
#
#for idx,row in blat.iterrows():
#    chr_ranges = exon_ranges['chr%s' % row.Tname]
#    
#    block_sizes = np.array([int(block) for block in row.blockSizes.split(',')[:-1]])
#    tstarts = np.array([int(start) for start in row.tStarts.split(',')[:-1]])
#    qstarts = np.array([int(start) for start in row.qStarts.split(',')[:-1]])
#    tends, qends = tstarts + block_sizes, qstarts + block_sizes
#    
#    for tstart, tend in zip(tstarts, tends):
#        
#        import ipdb; ipdb.set_trace()
