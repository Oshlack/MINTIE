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
from refine_annotations import get_pos_parts
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
    parser.add_argument(dest='vaf_estimates',
                        metavar='VAF_ESTIMATES',
                        type=str,
                        help='''VAF estimates file.''')
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
    contigs = contigs.drop_duplicates(ignore_index=True)

    if 'valid_motif' in contigs.columns:
        # remove NaNs (which are otherwise dropped when using groupby)
        contigs['valid_motif'] = contigs.valid_motif.fillna('Untested')

    return contigs

def get_short_gene_name(overlapping_genes):
    '''
    Extract first gene of list of overlapping
    genes at each fusion 'end'
    '''
    sgn = [og.split('|')[0] for og in overlapping_genes.split(':')]
    return '|'.join([g for g in sgn if g != ''])

def reformat_fields(contigs):
    '''
    Extract chrom, pos and strand fields.
    Reorder fields for clarity.
    Sort by p value.
    '''
    pos1 = contigs.pos1.apply(get_pos_parts).values
    pos2 = contigs.pos2.apply(get_pos_parts).values
    chr1, pos1, str1 = zip(*pos1)
    chr2, pos2, str2 = zip(*pos2)
    contigs['chr1'], contigs['pos1'], contigs['strand1'] = chr1, pos1, str1
    contigs['chr2'], contigs['pos2'], contigs['strand2'] = chr2, pos2, str2

    basic = ['chr1', 'pos1', 'strand1',
             'chr2', 'pos2', 'strand2',
             'variant_type', 'overlapping_genes', 'sample']
    variant = ['variant_id', 'partner_id', 'VAF', 'varsize',
               'contig_varsize', 'cpos', 'large_varsize',
               'is_contig_spliced', 'spliced_exon',
               'overlaps_exon', 'overlaps_gene']
    variant = variant if 'valid_motif' not in contigs.columns.values else variant + ['valid_motif']
    de = ['TPM', 'mean_WT_TPM']
    ec = ['ec_names', 'n_contigs_in_ec', 'contigs_in_EC', 'case_reads', 'controls_total_reads']
    cont = ['contig_id', 'unique_contig_ID', 'contig_len', 'contig_cigar']
    contigs = contigs[basic + variant + de + ec + cont]

    contigs = contigs.sort_values(by='case_reads', ascending=False)
    return contigs

def main():
    args = parse_args(sys.argv[1:])
    init_logging(args.log)

    try:
        contigs = pd.read_csv(args.contig_info, sep='\t', low_memory=False).fillna('')
        de_results = pd.read_csv(args.de_results, sep='\t', low_memory=False)
        vafs = pd.read_csv(args.vaf_estimates, sep='\t', low_memory=False)
        vafs = vafs[['contig_id', 'TPM', 'mean_WT_TPM', 'VAF']].drop_duplicates()
        gene_filter = pd.read_csv(args.gene_filter, header=None, low_memory=False) if args.gene_filter != '' \
                                                                                   else pd.DataFrame()
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    # consider only variants of interest
    contigs = contigs[contigs.variant_of_interest]
    contigs['sample'] = args.sample

    if args.var_filter:
        contigs = contigs[contigs.variant_type.apply(lambda v: v in args.var_filter).values]

    if len(gene_filter) > 0:
        contigs = filter_by_gene(contigs, gene_filter)

    if len(contigs) == 0:
        logging.info('WARNING: no variants present after filtering. Exiting.')
        contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')
        sys.exit()

    logging.info('Adding DE and VAF info...')
    contigs = add_de_info(contigs, de_results)
    contigs = pd.merge(contigs, vafs, on='contig_id', how='left')

    short_gnames = contigs.overlapping_genes.map(str).apply(get_short_gene_name)
    contig_ids, samples = contigs.contig_id, contigs['sample']
    con_names = ['|'.join([s, cid, sg]) for cid, s, sg in zip(contig_ids, samples, short_gnames)]
    contigs['unique_contig_ID'] = con_names

    logging.info('Outputting to CSV')
    contigs = reformat_fields(contigs)
    contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')

if __name__ == '__main__':
    main()
