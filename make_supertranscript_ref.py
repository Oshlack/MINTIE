#########################################################
# Author: Marek Cmero
#########################################################

import argparse
import pandas as pd
import numpy as np
import pysam
import os
import re
import gffutils
from Bio import SeqIO
from intervaltree import Interval, IntervalTree

import ipdb

parser = argparse.ArgumentParser()
parser.add_argument(dest='samfile',
                    help='''SAM or BAM format file containing contig alignments''')
parser.add_argument(dest='novel_contigs_file',
                    help='''Annotated novel contigs file''')
parser.add_argument(dest='blocks',
                    help='''Sequence blocks used to construct transcript sequences.''')
parser.add_argument(dest='fasta',
                    help='''Fasta file containing sequences (contructed from the blocks bed file)''')

args        = parser.parse_args()
samfile     = args.samfile
nc_file     = args.novel_contigs_file
blocks      = args.blocks
fasta_file  = args.fasta

novel_contigs = pd.read_csv(nc_file, sep='\t')
#novel_contigs = novel_contigs[novel_contigs.variant != 'soft-clip']
exon_df =  pd.read_csv(blocks, sep='\t', header = None, names = ['chrom', 'start', 'end', 'name', 'value', 'strand'])
exon_df[['gene', 'blocks']] = exon_df['name'].str.split(' ', expand=True)

seqs = {}
handle = open(fasta_file, 'r')
for record in SeqIO.parse(handle, 'fasta'):
    seqs[record.id] = str(record.seq)
handle.close()

bamf = pysam.AlignmentFile(samfile)
reads = pysam.IndexedReads(bamf)
reads.build()

gene_out = pd.DataFrame()
genes = np.unique([gn for gene in novel_contigs.gene.values for gn in gene.split('|')])
for gene in genes:
    nc = novel_contigs[novel_contigs.gene == gene]
    gene_out = exon_df[exon_df.gene == gene].copy()
    if len(gene_out) == 0:
        continue
    antisense = True if gene_out.strand.values[0]=='-' else False

    contigs_with_novel_bits = nc[nc.contig_varsize.values > 0]
    if any(nc.contig_varsize.values > 0):
        for idx, contig_row in contigs_with_novel_bits.iterrows():
            contig_mapping = reads.find(contig_row['contig'])
            seqs = [read.query_sequence for read in contig_mapping]
            #ipdb.set_trace()

    for idx,x in gene_out.iterrows():
        chrom = 'MT' if x.chrom == 'chrM' else x.chrom.split('chr')[1]
        seq = seqs['%s:%d-%d(%s)' % (chrom, x.start+1, x.end, x.strand)]

