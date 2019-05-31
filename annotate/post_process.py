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

# DE filters
MIN_LOGFC = 5

# score variables
VAR_WEIGHT = {'FUS': 1, 'INS': 1, 'DEL': 1, 'UN': 0.7, 'NE': 0.5, 'RI': 0.25, 'EE': 0.2, 'AS': 0.1}
EXP_WEIGHT = 0.6 # weight of expression related components; rest of weight is variant type

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
                        help='''Types of variant to keep.''')

    return parser.parse_args()

def get_all_genes(overlapping_genes):
    if isinstance(overlapping_genes, str):
        genes = overlapping_genes.split(':')
        genes = [gene.split('|') for gene in genes]
        genes = [g for gene in genes for g in gene if g != '']
        return genes
    else:
        return []

def filter_by_gene(contigs, gene_filter):
    genelist = gene_filter[0].values
    overlapping_genes = contigs.overlapping_genes.apply([lambda og: get_all_genes(og)])
    overlapping_genes = overlapping_genes.apply([lambda og: len(np.intersect1d(np.array(og), genelist)) > 0])
    contigs = contigs[overlapping_genes.values]
    return contigs

def add_de_info(contigs, de_results):
    de_results = de_results.rename(columns={'contig': 'contig_id', 'contigs': 'contigs_in_EC'})
    contigs = pd.merge(contigs, de_results, on='contig_id')
    contigs = contigs.drop(['genes'], axis=1)
    contigs = contigs[contigs.logFC >= MIN_LOGFC]

    # aggregate columns by variant ID
    agg_dict = {}
    for cname in contigs.columns.values:
        if cname in ['case_reads', 'controls_total_reads']:
            agg_dict[cname] = 'sum'
        elif cname in ['logFC', 'F', 'logCPM', 'n_txs_in_ec']:
            agg_dict[cname] = 'max'
        elif cname in ['contigs_in_EC', 'ec_names']:
            agg_dict[cname] = lambda x: ','.join(x)
        else:
            agg_dict[cname] = 'min'

    return contigs.groupby(by='contig_id').agg(agg_dict)

def get_st_alignments(contigs, st_bam):
    bam = pysam.AlignmentFile(st_bam, 'rc')
    st_alignment = []
    for contig in contigs.contig_id.values:
        aligned_conts = [read.reference_name for read in bam.fetch() if read.query_name == contig]
        aligned_conts = ','.join(np.unique(aligned_conts))
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
    #TODO: handle fusions
    contigs['crossing_reads'] = np.float('nan')
    contigs['junctions'] = np.float('nan')
    for idx,row in contigs.iterrows():
        st = '%s\|%s' % (row['sample'], row['contig_id'])
        st_blocks = st_bed[st_bed.contig.str.contains(st)]
        if len(st_blocks) > 0:
            st_blocks = make_junctions(st_blocks)
            rc = cjr.get_read_counts(read_align, st_blocks)
            contigs.loc[idx, 'crossing_reads'] = ','.join(rc.crossing.apply(str).values)
            se = zip(rc.start.values, rc.end.values)
            contigs.loc[idx, 'junctions'] = ','.join(['%s-%s' % (s,e) for s,e in se])
    return contigs

def add_score(contigs):
    '''
    add a score per variant based on:
    - expression metrics
        - p-value
        - logFC
        - reads in case EC(s)
        - reads in control EC(s)
    - variant type (based on pre-defined weights)
    '''
    top_pval = np.max(np.negative(np.log(contigs.PValue)))
    pval_score = np.negative(np.log(contigs.PValue)) / top_pval

    top_logFC = np.max(contigs.logFC)
    logFC_score = contigs.logFC / top_logFC

    case_read_score = contigs.case_reads / np.max(contigs.case_reads)
    con_read_score = 1 - (contigs.controls_total_reads / np.max(contigs.controls_total_reads))

    exp_score = 1/4 * np.array([pval_score, logFC_score, con_read_score, case_read_score])
    exp_score = np.sum(exp_score, axis=0) * EXP_WEIGHT

    tsv_score = np.array([VAR_WEIGHT[var] for var in contigs.variant_type]) * (1 - EXP_WEIGHT)

    scores = exp_score + tsv_score
    contigs['score'] = scores
    contigs = contigs.sort_values(by='score', ascending=False)

    return contigs

def main():
    args = parse_args()
    init_logging(args.log)

    try:
        contigs = pd.read_csv(args.contig_info, sep='\t', low_memory=False)
        de_results = pd.read_csv(args.de_results, sep='\t', low_memory=False)
        gene_filter = pd.read_csv(args.gene_filter, header=None, low_memory=False) if args.gene_filter != '' \
                                                                                   else pd.DataFrame()
        st_bed = pd.read_csv(args.st_bed, sep='\t', header=None, names=BED_COLS, low_memory=False)
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    # consider only variants of interest
    contigs = contigs[contigs.variant_of_interest]

    if args.var_filter:
        contigs = contigs[contigs.variant_type.apply(lambda v: v in args.var_filter).values]
        st_bed = st_bed[st_bed['name'].apply(lambda v: v in args.var_filter).values]

    contigs['sample'] = args.sample
    if len(gene_filter) > 0:
        contigs = filter_by_gene(contigs, gene_filter)

    contigs = add_de_info(contigs, de_results)
    contigs = get_st_alignments(contigs, args.cont_align)
    contigs = get_crossing_reads(contigs, args.read_align, st_bed)
    contigs = add_score(contigs)

    contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')

if __name__ == '__main__':
    main()
