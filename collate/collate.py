'''
Module      : collate
Description : Check ST alignment and count reads at coundaries
Copyright   : (c) Marek Cmero, Sep 2018
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
Features:
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
    parser.add_argument(dest='contig_info',
                        metavar='CONTIG_INFO',
                        type=str,
                        help='''Contig information for novel contigs.''')
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

    return parser.parse_args(args)

def get_short_gene_name(overlapping_genes):
    '''
    Extract first gene of list of overlapping
    genes at each fusion 'end'
    '''
    sgn = [og.split('|')[0] for og in overlapping_genes.split(':')]
    return '|'.join([g for g in sgn if g != ''])

def get_st_alignments(contigs, st_bam):
    bam = pysam.AlignmentFile(st_bam, 'rc')
    index = pysam.IndexedReads(bam)
    index.build()

    st_alignment = []
    for contig in contigs.contig_id.values:
        aligned_conts = [read.reference_name for read in index.find(contig) if read.reference_name]
        aligned_conts = ','.join(np.unique(aligned_conts)) if len(aligned_conts) > 0 else ''
        st_alignment.append(aligned_conts)

    # get short gene name (first gene of every overlapping set of genes, include fusion genes)
    short_gnames = contigs.overlapping_genes.apply(get_short_gene_name)
    contig_ids, samples = contigs.contig_id, contigs['sample']
    con_names = ['|'.join([s, cid, sg]) for cid, s, sg in zip(contig_ids, samples, short_gnames)]

    contigs['expected_ST_alignment'] = con_names
    contigs['real_ST_alignment'] = st_alignment
    return contigs

def make_junctions(st_blocks):
    for idx,row in st_blocks.iterrows():
        if row.end - row.start > SPLIT_LEN:
            start, end = row.start, row.end
            row['end'] = start
            st_blocks = st_blocks.append(row)
            row['start'], row['end'] = end, end
            st_blocks = st_blocks.append(row)
    return st_blocks[st_blocks.end - st_blocks.start <= SPLIT_LEN].drop_duplicates()

def get_read_support(contigs, bamf, st_bed):
    #TODO: handle fusions
    contigs['crossing_reads'] = np.float('nan')
    contigs['junctions'] = np.float('nan')
    for idx,row in contigs.iterrows():
        st = r'%s\|%s' % (row['sample'], row['contig_id'])
        st_blocks = st_bed[st_bed.contig.str.contains(st)]
        if len(st_blocks) > 0:
            st_blocks = make_junctions(st_blocks)
            rc = cjr.get_read_counts(bamf, st_blocks)
            contigs.loc[idx, 'crossing_reads'] = ','.join(rc.crossing.apply(str).values)
            se = zip(rc.start.values, rc.end.values)
            contigs.loc[idx, 'junctions'] = ','.join(['%s-%s' % (s,e) for s,e in se])
    return contigs

def main():
    args = parse_args(sys.argv[1:])
    init_logging(args.log)

    try:
        contigs = pd.read_csv(args.contig_info, sep='\t', low_memory=False)
        st_bed = pd.read_csv(args.st_bed, sep='\t', header=None, names=BED_COLS, low_memory=False)
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
    
    logging.info('Matching contigs to ST alignments...')
    contigs = get_st_alignments(contigs, args.cont_align)
    logging.info('Counting reads crossing variant boundaries...')
    bamf = pysam.AlignmentFile(args.read_align, "rb")
    contigs = get_read_support(contigs, bamf, st_bed)

    logging.info('Outputting to CSV')
    contigs = contigs.sort_values(by='PValue', ascending=True)
    contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')

if __name__ == '__main__':
    main()
