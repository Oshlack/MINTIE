#########################################################
# Author: Marek Cmero
# Take an samfile of aligned contigs, and filter out any
# contigs that match annotated transcript characteristics
# i.e. - genomic gaps <7bp
#      - no soft or hard clips >30bp
#      - all junctions are known (in annotation)
#########################################################

import argparse
import pandas as pd
import numpy as np
import pysam
import os

# cutoff parameters
gap_min = 7
clip_min = 30

# CIGAR specification codes
gaps = {'deletion': 1, 'insertion': 2, 'silent_deletion': 6}
clips = {'soft': 4, 'hard': 5}

def get_juncs(tx):
    '''
    return list of junctions in form
    [(chr, start, end)] from transcript
    info file
    '''
    starts = tx['exonStarts'].split(',')[1:-1]
    ends = tx['exonEnds'].split(',')[:-2]
    chroms = [tx['chrom']] * len(starts)
    return(list(zip(chroms, ends, starts)))

parser = argparse.ArgumentParser()
parser.add_argument(dest='tx_info')
parser.add_argument(dest='samfile')
parser.add_argument(dest='outbam_file')
args = parser.parse_args()

tx_info = args.tx_info
samfile = args.samfile
outbam_file = args.outbam_file

genref = pd.read_csv(tx_info, sep='\t')
juncs = genref.apply(lambda tx: get_juncs(tx), axis=1)
juncs = [(str(c), int(s), int(e)) for jv in juncs.values for c, s, e in jv] # flatten juncs list

sam = pysam.AlignmentFile(samfile, 'rc')
outbam = pysam.AlignmentFile(outbam_file, 'wb', template=sam)

# write novel contigs to bam file
int_contigs = np.empty(0, dtype='U100')
for read in sam.fetch():
    has_gaps = any([op in gaps.values() and val >= gap_min for op, val in read.cigar])
    has_clips = any([op in clips.values() and val >= clip_min for op, val in read.cigar])

    # check junctions
    starts, ends = zip(*read.blocks)
    chroms = [read.reference_name] * (len(starts)-1)
    tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
    unknown_juncs = any([txj not in juncs for txj in tx_juncs])

    if has_gaps or has_clips or unknown_juncs:
        int_contigs = np.append(int_contigs, read.qname)
        outbam.write(read)

# write interesting contigs list to file
int_contigs = np.unique(int_contigs)
outdir = os.path.dirname(outbam_file)
outdir = '.' if outdir == '' else outdir
with open('%s/interesting_contigs.txt' % outdir, 'w') as fout:
    for contig in int_contigs:
        fout.write('%s\n' % contig)
