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

def split_block(blocks, block, block_seqs, gpos1, gpos2, seq, name, strand):
    '''
    Split a sequence block at the gpos location, separating
    into left and right blocks and sequences (if necessary)
    '''
    blocks = blocks[blocks.index!=block['index']]
    ref_seq = block_seqs['%s:%d-%d(%s)' % (block['chr'], block.start, block.end, strand)]
    end = -1 if 'DEL' else len(ref_seq) - 1
    left_seq = ref_seq[1:gpos1-block.start]
    right_seq = ref_seq[gpos2-block.start-1:end]

    block_seqs['%s:%d-%d(%s)' % (block['chr'], gpos1, gpos2, strand)] = str(seq)
    var_block = [{'chr': block['chr'], 'start': gpos1, 'end': gpos2, \
                  'name': name, 'score': '.', 'strand': strand}]
    blocks = blocks.append(var_block, ignore_index=True)

    if gpos1 - block.start > 0:
        block_seqs['%s:%d-%d(%s)' % (block['chr'], block.start, gpos1, strand)] = left_seq
        left_block = [{'chr': block['chr'], 'start': block.start, 'end': gpos1, \
                       'name': block['name'], 'score': '.', 'strand': strand}]
        blocks = blocks.append(left_block, ignore_index=True)

    if block.end - gpos2 > 0:
        block_seqs['%s:%d-%d(%s)' % (block['chr'], gpos2, block.end, strand)] = right_seq
        right_block = [{'chr': block['chr'], 'start': gpos2, 'end': block.end, \
                        'name': block['name'], 'score': '.', 'strand': strand}]
        blocks = blocks.append(right_block, ignore_index=True)

    return blocks, block_seqs

def sort_blocks(blocks):
    '''
    Sort blocks in ascending order if on the sense strand,
    and descending order if on the antisense strand
    '''
    blocks = blocks.drop_duplicates().sort_values(by=['chr','start','end']).reset_index(drop=True)
    if '-' in blocks.strand.values:
        antisense_blocks = blocks[blocks.strand=='-'].sort_values(by=['chr','start','end'], ascending=False)
        sense_blocks = blocks[blocks.strand=='+']
        if len(sense_blocks) > 0:
            sense_first = blocks.strand.values[0] == '+'
            blocks = sense_blocks.append(antisense_blocks) \
                        if sense_first \
                        else antisense_blocks.append(sense_blocks)
        else:
            blocks = antisense_blocks
    return blocks

def get_block_colours(blocks, names, alt=False):
    '''
    Get alternating colours for supertranscript block bed
    '''
    colours = np.empty((len(blocks),), dtype='U50')
    colours[::2] = COL1 if not alt else COL3
    colours[1::2] = COL2 if not alt else COL4
    novel_vars = [x in VARS_TO_ANNOTATE for x in names]
    colours[novel_vars] = VARCOL
    return colours

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
    if len(left) > 0 and len(right) > 0 and block_size > MIN_CLIP:
        return True

    if len(left) > 0:
        block_size = left.start.values[0] - block[0]
    elif len(right) > 0:
        block_size = block[1] - right.end.values[0]

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
