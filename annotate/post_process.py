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
from Bio import SeqIO
from argparse import ArgumentParser
from utils import cached, init_logging, exit_with_error
import ipdb

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1

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
    #parser.add_argument(dest='st_align_bam',
    #                    metavar='ST_ALIGN_BAM',
    #                    type=str,
    #                    help='''Alignment of novel contigs to supertranscript reference.''')
    parser.add_argument('--gene_filter',
                        metavar='GENE_FILTER',
                        type=str,
                        default='',
                        help='''List of genes (one per line) to keep (filter out others).''')

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
    return contigs

def main():
    args = parse_args()
    init_logging(args.log)

    try:
        contigs = pd.read_csv(args.contig_info, sep='\t')
        de_results = pd.read_csv(args.de_results, sep='\t')
        gene_filter = pd.read_csv(args.gene_filter, header=None) if args.gene_filter != '' else pd.DataFrame()
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    contigs['sample'] = args.sample
    if len(gene_filter) > 0:
        contigs = filter_by_gene(contigs, gene_filter)

    contigs = add_de_info(contigs, de_results)

if __name__ == '__main__':
    main()
