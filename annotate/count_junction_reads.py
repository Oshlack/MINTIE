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
import logging

rc_dtype = [('contig', '<U150'),
            ('start', 'int64'),
            ('end', 'int64'),
            ('crossing', 'int64')]

def get_crossing_reads(contig, start, end, bamf):
    try:
        crossing_reads = []
        iter_loc = bamf.fetch(contig, start, end, until_eof=True)
        for read in iter_loc:
            if read.reference_start < start and read.reference_end > end:
                crossing_reads.append(read.query_name)
        return np.unique(crossing_reads)
    except ValueError:
        loc = '%s:%d-%d' % (contig, start, end)
        logging.info('Fetching reads failed at loc %s; skipping.' % loc)
        return np.array([])

def get_read_counts(bam_path, juncs):
    bamf = pysam.AlignmentFile(bam_path, "rb")
    read_counts = np.empty([0, len(rc_dtype)], dtype=rc_dtype)
    for idx, junc in juncs.iterrows():
        contig = junc['contig']
        start = junc['start']
        end = junc['end']+1

        loc = '%s:%d-%d' % (contig, start, end)
        crossing_reads = len(get_crossing_reads(contig, start, end, bamf))

        junc_rc = np.array([(contig, start, end, crossing_reads)], dtype=rc_dtype)
        read_counts = np.append(read_counts, junc_rc)

    return pd.DataFrame(read_counts)
