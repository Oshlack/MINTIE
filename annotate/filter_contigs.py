
'''
Module      : filter_contigs
Description : Filter out any contigs matching reference
Copyright   : (c) Marek Cmero, Sep 2018
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
Take an samfile of aligned contigs and filter out
any contigs matching the reference annotation.
'''

import pandas as pd
import numpy as np
import pysam
import os
import logging
import ipdb
import annotate_contigs as ac
from argparse import ArgumentParser
from utils import cached, init_logging, exit_with_error

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "annotate_contigs"

# CUT-OFF parameters
GAP_MIN = 7
CLIP_MIN = 30
MATCH_MIN = 30
MATCH_PERC_MIN = 0.3

# VCF parameters
INFO = ["CID", "ECN", "CLEN", "CPOS", "CSTRAND", "CCIGAR", "VSIZE", "CVSIZE", "CVTYPE", "GENES", "PARID", "PVAL", "CVQ"]
FORMAT = ["GT", "ECC", "AI"]

# CIGAR specification codes
CIGAR = {'match': 0,
         'insertion': 1,
         'deletion': 2,
         'skipped': 3,
         'soft-clip': 4,
         'hard-clip': 5,
         'silent_deletion': 6}

GAPS = [CIGAR[c] for c in ['insertion', 'deletion', 'silent_deletion']]
CLIPS = [CIGAR[c] for c in ['soft-clip', 'hard-clip']]
# any cigar criteria that is >0 bp on an aligned contig
AFFECT_CONTIG = [CIGAR[c] for c in ['match', 'insertion', 'soft-clip', 'hard-clip']]
# any cigar criteria that is >0 bp on the reference genome
AFFECT_REF = [CIGAR[c] for c in ['match', 'deletion']]

# read processing record
record = {}

def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Filter contigs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='sam_file',
                        metavar='SAM_FILE',
                        type=str,
                        help='''SAM or BAM format file containing contig alignments''')
    parser.add_argument(dest='junc_file',
                        metavar='JUNC_FILE',
                        type=str,
                        help='''Reference file containing transcripts and their
                        respective splice junctions.''')
    parser.add_argument(dest='tx_ref_file',
                        type=str,
                        metavar='TX_REF_FILE',
                        help='''Transcriptiome GTF reference file.''')
    return parser.parse_args()

def filter_contigs(args):
    '''
    Iterate through contig reads of bam file and
    outout the names of any contigs that do not
    match the reference
    '''
    sam = pysam.AlignmentFile(args.sam_file, 'rc')
    ref_trees, ex_ref = ac.get_gene_lookup(args.tx_ref_file)
    juncs = ac.get_junc_lookup(args.junc_file)

    logging.info('Checking contigs for non-reference content...')
    novel_contigs = []
    for read in sam.fetch():
        genes = ac.get_overlapping_genes(read, ref_trees)
        if genes == '':
            logging.info('No gene(s) intersecting read %s; skipping' % read.query_name)
            continue

        # check for contig gaps or soft-clips
        has_gaps = any([op in GAPS and val >= GAP_MIN for op, val in read.cigar])
        has_clips = any([op in CLIPS and val >= CLIP_MIN for op, val in read.cigar])

        # check junctions
        tx_juncs = ac.get_tx_juncs(read)
        unknown_juncs = ['%s:%s-%s' % (c, s, e) not in juncs[0] for c, s, e in tx_juncs]
        has_novel_juncs = any(unknown_juncs)

        # check for novel blocks
        chr_ref = ac.get_chr_ref(read, ex_ref)
        has_novel_blocks = any([ac.is_novel_block(block, chr_ref) for block in read.get_blocks()])

        if has_gaps or has_clips or has_novel_juncs or has_novel_blocks:
            print("\t".join([read.query_name, genes]))

def main():
    args = parse_args()
    init_logging(args.log)
    filter_contigs(args)

if __name__ == '__main__':
    main()
