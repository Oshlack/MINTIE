'''
Module      : block_helper
Description : Helper functions for handling exon blocks
Copyright   : (c) Marek Cmero, Mar 2019
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
'''
import numpy as np

VARS_TO_ANNOTATE = ['EE','NE','INS','RI','UN','FUS']

# alternating colours for bed track, and variant colour
COL1 = '99,99,99' # dark grey
COL2 = '189,189,189' #light grey
COL3 = '49,130,189' # dark blue
COL4 = '222,235,247' # light blue
VARCOL = '255,255,153' # bright yellow

def split_block(blocks, block, block_seqs, gpos1, gpos2, seq, name, strand):
    '''
    split a sequence block at the gpos location, separating
    into left and right blocks and sequences (if necessary)
    '''
    blocks = blocks[blocks.index!=block['index']]
    ref_seq = block_seqs['%s:%d-%d(%s)' % (block['chr'], block.start, block.end, strand)]
    left_seq = ref_seq[:gpos1-block.start]
    right_seq = ref_seq[gpos2-block.start-1:]

    block_seqs['%s:%d-%d(%s)' % (block['chr'], gpos1, gpos2, strand)] = str(seq)
    var_block = [{'chr': block['chr'], 'start': gpos1, 'end': gpos2, \
                  'name': name, 'score': '.', 'strand': strand}]
    blocks = blocks.append(var_block, ignore_index=True)

    if gpos1 - block.start != 0:
        block_seqs['%s:%d-%d(%s)' % (block['chr'], block.start, gpos1, strand)] = left_seq
        left_block = [{'chr': block['chr'], 'start': block.start, 'end': gpos1, \
                       'name': block['name'], 'score': '.', 'strand': strand}]
        blocks = blocks.append(left_block, ignore_index=True)

    if block.end - gpos2 != 0:
        block_seqs['%s:%d-%d(%s)' % (block['chr'], gpos2, block.end, strand)] = right_seq
        right_block = [{'chr': block['chr'], 'start': gpos2, 'end': block.end, \
                        'name': block['name'], 'score': '.', 'strand': strand}]
        blocks = blocks.append(right_block, ignore_index=True)

    return blocks, block_seqs

def sort_blocks(blocks):
    '''
    sort blocks in ascending order if on the sense strand,
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
    colours = np.empty((len(blocks),), dtype='U50')
    colours[::2] = COL1 if not alt else COL3
    colours[1::2] = COL2 if not alt else COL4
    novel_vars = [x in VARS_TO_ANNOTATE for x in names]
    colours[novel_vars] = VARCOL
    return colours

