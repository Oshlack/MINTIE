'''
Module      : block_helper
Description : Helper functions for handling exon blocks
Copyright   : (c) Marek Cmero, Mar 2019
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
'''
import numpy as np
import constants

VARS_TO_ANNOTATE = ['EE','NE','INS','RI','UN','FUS','DEL']

# alternating colours for bed track, and variant colour
COL1 = '99,99,99' # dark grey
COL2 = '189,189,189' #light grey
COL3 = '49,130,189' # dark blue
COL4 = '222,235,247' # light blue
VARCOL = '255,255,153' # bright yellow

def is_novel_block(block, chr_ref, MIN_CLIP):
    '''
    Checks a contig's sequence blocks and returns
    false if the block matches (or is contained)
    within a referenced exon, and true otherwise
    '''
    contained = chr_ref[np.logical_and(block[0] >= chr_ref.start, block[1] <= chr_ref.end)]
    if len(contained) > 0:
        return False

    block_size = block[1] - block[0]
    olapping = chr_ref[np.logical_and(block[0] < chr_ref.start, block[1] > chr_ref.end)]
    if len(olapping) > 0 and block_size > MIN_CLIP:
        return True

    left = chr_ref[np.logical_and(block[1] > chr_ref.start, block[1] <= chr_ref.end)]
    right = chr_ref[np.logical_and(block[0] >= chr_ref.start, block[0] < chr_ref.end)]
    if len(left) > 0 and len(right) > 0:
        block_size = min(left.start.values) - max(right.end.values)
        assert block_size >= 0
        return block_size >= MIN_CLIP

    if len(left) > 0:
        block_size = min(left.start.values) - block[0]
    elif len(right) > 0:
        block_size = block[1] - max(right.end.values)

    assert block_size >= 0
    return block_size >= MIN_CLIP

def get_block_sequence(read, block_idx):
    '''
    Return query and reference sequence of
    specified block from contig read
    '''
    qseq = read.query_sequence
    rseq = read.get_reference_sequence()

    cpos1 = sum([v for c,v in read.cigar[:block_idx] if c in constants.AFFECT_CONTIG and c != constants.CIGAR['hard-clip']])
    cpos2 = sum([v for c,v in read.cigar[:block_idx+1] if c in constants.AFFECT_CONTIG and c != constants.CIGAR['hard-clip']])

    rpos1 = sum([v for c,v in read.cigar[:block_idx] if c in constants.AFFECT_REF])
    rpos2 = sum([v for c,v in read.cigar[:block_idx+1] if c in constants.AFFECT_REF])

    return qseq[cpos1:cpos2], rseq[rpos1:rpos2]
