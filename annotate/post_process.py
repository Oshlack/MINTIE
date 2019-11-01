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
from argparse import ArgumentParser
from utils import init_logging, exit_with_error

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1
BED_COLS = ['contig', 'start', 'end', 'name', 'score', 'strand', 'tStart', 'tEnd', 'itemRgb']
SPLIT_LEN = 10 # split variants longer than this many base-pairs into two separate junctions to count reads for

def parse_args(args):
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

    return parser.parse_args(args)

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

    # aggregate columns by variant ID
    agg_dict = {}
    for cname in contigs.columns.values:
        if cname in ['case_reads', 'controls_total_reads']:
            agg_dict[cname] = 'sum'
        elif cname in ['logFC', 'F', 'logCPM', 'n_contigs_in_ec']:
            agg_dict[cname] = 'max'
        elif cname in ['PValue', 'FDR']:
            agg_dict[cname] = 'min'
        elif cname in ['contigs_in_EC', 'ec_names']:
            agg_dict[cname] = lambda x: ','.join(x)

    group_vars = [c for c in contigs.columns if c not in agg_dict.keys()]
    return contigs.groupby(by=group_vars, as_index=False).agg(agg_dict)

def get_short_gene_name(overlapping_genes):
    '''
    Extract first gene of list of overlapping
    genes at each fusion 'end'
    '''
    sgn = [og.split('|')[0] for og in overlapping_genes.split(':')]
    return '|'.join([g for g in sgn if g != ''])

def main():
    args = parse_args(sys.argv[1:])
    init_logging(args.log)

    try:
        contigs = pd.read_csv(args.contig_info, sep='\t', low_memory=False)
        de_results = pd.read_csv(args.de_results, sep='\t', low_memory=False)
        gene_filter = pd.read_csv(args.gene_filter, header=None, low_memory=False) if args.gene_filter != '' \
                                                                                   else pd.DataFrame()
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    # consider only variants of interest
    contigs = contigs[contigs.variant_of_interest]

    if args.var_filter:
        contigs = contigs[contigs.variant_type.apply(lambda v: v in args.var_filter).values]

    contigs['sample'] = args.sample
    if len(gene_filter) > 0:
        contigs = filter_by_gene(contigs, gene_filter)

    logging.info('Adding DE info...')
    contigs = add_de_info(contigs, de_results)

    short_gnames = contigs.overlapping_genes.apply(get_short_gene_name)
    contig_ids, samples = contigs.contig_id, contigs['sample']
    con_names = ['|'.join([s, cid, sg]) for cid, s, sg in zip(contig_ids, samples, short_gnames)]
    contigs['unique_contig_ID'] = con_names

    logging.info('Outputting to CSV')
    contigs = contigs.sort_values(by='PValue', ascending=True)
    contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')

if __name__ == '__main__':
    main()
