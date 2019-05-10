'''
Module      : post_process
Description : Filter and collate novel variant info
Copyright   : (c) Marek Cmero, Sep 2018
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
Features:
    - collate expression and annotation information
    - filter by variant type
    - allow gene list filtering
    - check contig alignment to supertranscript
    - check reads crossing boundary of novel variant
'''
import numpy as np
import pandas as pd
import re
import sys
import logging
import os
import pysam
import count_junction_reads as cjr
from Bio import SeqIO
from argparse import ArgumentParser
from utils import cached, init_logging, exit_with_error
import ipdb

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1
BED_COLS = ['contig', 'start', 'end', 'name', 'score', 'strand', 'tStart', 'tEnd', 'itemRgb']
SPLIT_LEN = 10 # split variants longer than this many base-pairs into two separate junctions to count reads for

def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Make supertranscript reference'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='sample',
                        metavar='SAMPLE',
                        type=str,
                        help='''Sample name.''')
    parser.add_argument(dest='contig_info',
                        metavar='CONTIG_INFO',
                        type=str,
                        help='''Contig information for novel contigs.''')
    parser.add_argument(dest='de_results',
                        metavar='DE_RESULTS',
                        type=str,
                        help='''Differential expression results.''')
    parser.add_argument(dest='st_bed',
                        metavar='ST_BED',
                        type=str,
                        help='''Supertranscript variant bed file.''')
    parser.add_argument(dest='cont_align',
                        metavar='CONT_ALIGN',
                        type=str,
                        help='''Alignment of novel contigs to supertranscript reference.''')
    parser.add_argument(dest='read_align',
                        metavar='READ_ALIGN',
                        type=str,
                        help='''Alignment of reads to supertranscript reference.''')
    parser.add_argument('--gene_filter',
                        metavar='GENE_FILTER',
                        type=str,
                        default='',
                        help='''File containing list of genes (one per line) to keep (filter out others).''')
    parser.add_argument('--var_filter',
                        metavar='VAR_FILTER',
                        type=str,
                        nargs='+',
                        help='''Variants to keep.''')

    return parser.parse_args()

def get_all_genes(overlapping_genes):
    genes = overlapping_genes.split(':')
    genes = [gene.split('|') for gene in genes]
    genes = [g for gene in genes for g in gene]
    return genes

def is_in_genelist(overlapping_genes, genelist):
    return any([gene in genelist for gene in overlapping_genes])

def filter_by_gene(contigs, gene_filter):
    genelist = gene_filter[0].values
    overlapping_genes = contigs.overlapping_genes.apply([lambda og: get_all_genes(og)])
    overlapping_genes = overlapping_genes.apply([lambda og: is_in_genelist(og, genelist)])
    contigs = contigs[overlapping_genes.values]
    return contigs

def add_de_info(contigs, de_results):
    de_results = de_results.rename(columns={'contig': 'contig_id', 'contigs': 'contigs_in_EC'})
    contigs = pd.merge(contigs, de_results, on='contig_id')
    contigs = contigs.drop(['genes'], axis=1)
    return contigs

def get_st_alignments(contigs, st_bam):
    bam = pysam.AlignmentFile(st_bam, 'rc')
    st_alignment = []
    for contig in contigs.contig_id.values:
        aligned_conts = [read.reference_name for read in bam.fetch() if read.query_name == contig]
        aligned_conts = ','.join(aligned_conts)
        st_alignment.append(aligned_conts)
    contigs['ST_alignment'] = st_alignment
    return contigs

def make_junctions(st_blocks):
    for idx,row in st_blocks.iterrows():
        if row.end - row.start > SPLIT_LEN:
            start, end = row.start, row.end
            row['end'] = start
            st_blocks = st_blocks.append(row)
            row['start'], row['end'] = end, end
            st_blocks = st_blocks.append(row)
    return st_blocks[st_blocks.end - st_blocks.start <= SPLIT_LEN]

def get_crossing_reads(contigs, read_align, st_bed):
    #TODO: handle fusions and deletions
    contigs['crossing_reads'] = np.float('nan')
    contigs['junctions'] = np.float('nan')
    for idx,row in contigs.iterrows():
        st = row['ST_alignment']
        st_blocks = st_bed[np.logical_and(st_bed.contig == st, st_bed['name'] == row.variant_type)]
        if len(st_blocks) > 0:
            st_blocks = make_junctions(st_blocks)
            rc = cjr.get_read_counts(read_align, st_blocks)
            contigs.loc[idx, 'crossing_reads'] = ','.join(rc.crossing.apply(str).values)
            contigs.loc[idx, 'junctions'] = ','.join(['%s-%s' % (s,e) for s,e in zip(rc.start.values, rc.end.values)])
    return contigs

def main():
    args = parse_args()
    init_logging(args.log)

    try:
        contigs = pd.read_csv(args.contig_info, sep='\t')
        de_results = pd.read_csv(args.de_results, sep='\t')
        gene_filter = pd.read_csv(args.gene_filter, header=None) if args.gene_filter != '' else pd.DataFrame()
        st_bed = pd.read_csv(args.st_bed, sep='\t', header=None, names=BED_COLS)
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    if args.var_filter:
        contigs = contigs[contigs.variant_type.apply(lambda v: v in args.var_filter).values]
        st_bed = st_bed[st_bed['name'].apply(lambda v: v in args.var_filter).values]

    contigs['sample'] = args.sample
    if len(gene_filter) > 0:
        contigs = filter_by_gene(contigs, gene_filter)

    contigs = add_de_info(contigs, de_results)
    contigs = get_st_alignments(contigs, args.cont_align)
    contigs = get_crossing_reads(contigs, args.read_align, st_bed)

    contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')

if __name__ == '__main__':
    main()
