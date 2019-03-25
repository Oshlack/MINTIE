########################################################
# Author: Marek Cmero
# Count reads crossing block boundaries
########################################################

import argparse
import sys
import pandas as pd
import numpy as np
import pysam
import os
import re
import ipdb

read_dtype =   [('query_name', '<U150'),
                ('chrom', '<U50'),
                ('ref_start', 'int64'),
                ('ref_end', 'int64'),
                ('align_start', 'int64'),
                ('align_end', 'int64'),
                ('len', 'int64'),
                ('ins_len', 'int64'),
                ('is_reverse', np.bool)]

rc_dtype = [('contig', '<U150'),
            ('start', 'int64'),
            ('end', 'int64'),
            ('crossing', 'int64'),
            ('depth', 'int64')]

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_to_array(x,bamf):
    chrom = bamf.getrname(x.reference_id)
    try:
        read = np.array((x.query_name,chrom, x.reference_start, x.reference_end,
                         x.query_alignment_start, x.query_alignment_end,
                         x.query_length, x.tlen, np.bool(x.is_reverse)),
                        dtype=read_dtype)
        return read
    except TypeError:
        eprint('Warning: record %s contains invalid attributes, skipping' % x.query_name)
        return np.empty(0)

def get_loc_reads(contig, start, end, bamf, max_dp=10000):
    loc = '%s:%d-%d' % (contig, start, end)
    loc_reads = np.empty([0, len(read_dtype)], dtype=read_dtype)
    err_code = 0
    try:
        iter_loc = bamf.fetch(contig, start, end, until_eof=True)
        for x in iter_loc:
            read = read_to_array(x, bamf)
            if len(np.atleast_1d(read))>0:
                loc_reads = np.append(loc_reads, read)
            if len(loc_reads) > max_dp:
                eprint('Read depth too high at %s' % loc)
                err_code = 1
                return np.empty(0), err_code
        loc_reads = np.sort(loc_reads, axis=0, order=['query_name','ref_start'])
        loc_reads = np.unique(loc_reads) #remove duplicates
        return loc_reads, err_code
    except ValueError:        
        eprint('Fetching reads failed for loc: %s' % loc)
        err_code = 2
        return loc_reads, err_code

def read_pair_crosses_junction(read, mate, start, end):
    if mate == None:
        return (read['ref_start'] < start and read['ref_end'] > end)
    else:
        return (read['ref_start'] < start and mate['ref_end'] > end)
    
def get_read_counts(bam_path, juncs):
    bamf = pysam.AlignmentFile(bam_path, "rb")
    read_counts = np.empty([0, 5], dtype=rc_dtype)
    for idx, junc in juncs.iterrows():
        contig = junc['contig']
        start = junc['start']
        end = junc['end']+1

        loc = '%s:%d-%d' % (contig, start, end)
        loc_reads, err_code = get_loc_reads(contig, start, end, bamf)

        nreads = len(loc_reads)
        crossing_reads = 0

        for ridx, read in enumerate(loc_reads):
            r1 = loc_reads[ridx]
            r2 = loc_reads[ridx+1] if (ridx+2)<=len(loc_reads) else None
            if read_pair_crosses_junction(r1, r2, start, end):
                crossing_reads += 1

        junc_rc = np.array([(contig, start, end, crossing_reads, nreads)], dtype=rc_dtype)
        read_counts = np.append(read_counts, junc_rc)

    return pd.DataFrame(read_counts)
