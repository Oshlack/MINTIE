'''
Module      : simu
Description : Functions to help in running Jupyter notebook simulations
Copyright   : (c) Marek Cmero, Mar 2019
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
'''
import tempfile
import numpy as np
import pandas as pd
import collections
import re
import pybedtools
from IPython.core.debugger import set_trace
from intervaltree import Interval, IntervalTree
from pybedtools import BedTool
from Bio import SeqIO

BASES = list('GCAT')
SEED_INIT = 123

def get_gene_name(row):
    '''
    Prevents KeyError if gene name missing
    '''
    try:
        return row.attrs['gene_name']
    except KeyError:
        return ''
    
def get_seq(gr, genome_fasta):
    '''
    Returns a dictionary of exon sequences
    and the corresponding strand of the transcript
    '''
    block_seqs = gr.sequence(fi=genome_fasta, s=True)
    block_dict = collections.OrderedDict()
    with tempfile.NamedTemporaryFile() as fa_tmp:
        fa_tmp.write(bytes(open(block_seqs.seqfn).read(), 'utf-8'))
        fa_tmp.flush()

        for record in SeqIO.parse(fa_tmp.name, 'fasta'):
            block_dict[record.id] = str(record.seq)

    strand = re.search('\(([-+])\)', next(iter(block_dict.keys())))
    assert strand
    
    return block_dict, strand.group(1)

def write_sequence(seq_dict, strand, output_file, name):
    '''
    Writes sequence dictionary to output file
    '''
    seq = [seq_dict[ex] for ex in seq_dict.keys()]
    seq = seq if strand == '+' else [s for s in reversed(seq)]

    with open(output_file, 'a') as fout:
        fout.write('>%s\n' % name)
        fout.write(''.join(seq) + '\n')
        
def get_chrom_features(chrom, gr):
    '''
    Get all merged intervals on given chromosome
    '''
    chrom_features = gr.filter(lambda x: x.chrom == chrom).merge()    
    chrom_features = [(g.start, g.end) for g in gr]
    
    chrom_tree = IntervalTree()
    [chrom_tree.addi(s, e) for s, e in chrom_features]

    return chrom_tree

def get_gene_features(gr):
    '''
    Get interval tree for start/ends for 
    each gene on each chromosome
    '''
    chroms = np.unique([x.chrom for x in gr])
    gn_ref = pd.DataFrame([(g.chrom, g.start, g.end, get_gene_name(g)) for g in gr])
    aggregator = {1: lambda x: min(x),
                  2: lambda x: max(x)}
    gn_ref = gn_ref.groupby([0, 3], as_index=False, sort=False).agg(aggregator)
    gn_ref = gn_ref[[0, 1, 2, 3]]
    gn_ref.columns = ['chrom', 'start', 'end', 'gene']
    gn_ref = gn_ref[gn_ref.gene!='']

    ref_trees = {}
    for chrom in chroms:
        chr_ref = gn_ref[gn_ref.chrom == chrom]
        ref_tree = IntervalTree()
        for s,e,g in zip(chr_ref['start'].values, chr_ref['end'].values, chr_ref['gene'].values):
            ref_tree.addi(s-1, e, g)
    ref_trees[chrom] = ref_tree
    
    return ref_trees

def get_exon_seq(ex_list, strand, gr, genome_fasta, block_range, extended=True):
    '''
    Extends given exon, or creates a novel downstream
    exon with random size and returns its sequence. 
    If a reference exon exists that already extends
    or overlaps the given exon, it will extend or place
    the exon past the overlapping exon.
    '''
    # TODO: make sure exon size doesn't extend to the next exon 
    block_size = np.random.randint(block_range[0], block_range[1])
    gap_size = 0 if extended else np.random.randint(block_range[0], block_range[1])
    
    block = ex_list[1] if strand == '+' else ex_list[0]
    loc = re.compile('[:\-\(\)]').split(block)
    exon_start, exon_end = int(loc[1]), int(loc[2])
        
    start = exon_end + gap_size if strand == '+' else exon_start - gap_size - block_size
    end = start + block_size

    ex = get_chrom_features(loc[0], gr.merge())
    olap = ex.overlap(int(start), int(end))
    if len(olap) > 0:
        coords = list(olap)[0]
        s, e = int(coords[0]), int(coords[1])        
        if extended:
            start = start if strand == '+' else s - block_size
            end = end if strand == '-' else e + block_size
        else:          
            start = e + gap_size if strand == '+' else s - gap_size - block_size            
            end = start + block_size
    
    block_bed = '%s\t%d\t%d\t.\t1\t%s' % (loc[0], start, end, strand)
    block_bt = BedTool(block_bed, from_string=True)
    block_seq, bs = get_seq(block_bt, genome_fasta)
    ext_seq = ''.join([bs for bs in block_seq.values()])
    bloc = ''.join([k for k in block_seq.keys()])
    
    return ext_seq, bloc

def increment_seed(seed, amount=1):
    seed += amount
    np.random.seed(seed)
    return seed

def get_random_block(all_exons, gene_trees, genome_fasta, block_range):
    '''
    Get random block sequence for feature
    not overlapping any other genomic features
    '''
    chr_sizes = pybedtools.chromsizes('hg38')
    chroms = np.unique([x.chrom for x in all_exons])
    
    block_size = np.random.randint(block_range[0], block_range[1])
    chrom = np.random.choice(chroms)
    chrom_features = gene_trees[chrom]
    
    chr_range = chr_sizes[('chr%s' % chrom)]
    block_start = np.random.randint(chr_range[0], chr_range[1]-block_size)
    block_end = block_start + block_size        
    
    seed = SEED_INIT
    seq = 'N'
    while 'N' in seq:
        # only select sequence if there's no Ns
        while chrom_features.overlaps(block_start, block_end):
            block_start = np.random.randint(chr_range[0], chr_range[1]-block_size)
            block_end = block_start + block_size
            
            if chrom_features.overlaps(block_start, block_end):
                seed = increment_seed(seed)

        strand = np.random.choice(['+','-'])
        block_bed = '%s\t%d\t%d\t.\t1\t%s' % (chrom, block_start, block_end, strand)
        block_bt = BedTool(block_bed, from_string=True)
        block_seq, bs = get_seq(block_bt, genome_fasta)
        seq = ''.join([bs for bs in block_seq.values()])
        
        if 'N' in seq:
            seed = increment_seed(seed)

    return block_seq, seed

def get_random_seq(ins_range):
    '''
    Generate random insertion sequence
    '''
    ins_size = np.random.randint(ins_range[0], ins_range[1])
    ins = np.random.choice(BASES, ins_size)
    ins = ''.join(ins)
    return ins

def get_tx_seq(tx, all_exons, genome_fasta, n_exons, control_fasta, front=True, wt_out=True):
    '''
    Get fusion sequence of given transcript, returning
    sequence of first N exons for transcript 1 (front=True)
    and N exons for transcript 2 (front=False). By default,
    write out wildtype transcript to control reference file.
    '''
    
    exons = all_exons.filter(lambda x: x['transcript_id'] == tx).saveas()
    tx_seq, s = get_seq(exons, genome_fasta)
    ex_list = [ex for ex in tx_seq.keys()]
    
    # pick N 5' exons for tx1 (front) and N 3' exons for tx2 (back)
    if front:
        ex_list = ex_list[:n_exons] if s == '+' else ex_list[-n_exons:]
    else:
        ex_list = ex_list[-n_exons:] if s == '+' else ex_list[:n_exons]
    
    # select sequences and reverse order if antisense
    seq = [tx_seq[ex] for ex in ex_list]
    seq = seq if s == '+' else [s for s in reversed(seq)]  
    
    if wt_out:
        write_sequence(tx_seq, s, control_fasta, tx)
    
    return seq, s, ex_list

def write_fusion(tx1, tx2, all_exons, genome_fasta, params, gene_trees, add=None):
    '''
    Get left and right sequences of given transcripts
    corresponding to the first N exons and last N exons of
    transcripts 1 and 2 respectively.
    Automatically writes wild type transcript.
    '''    
    exon_types = ['EE', 'NE', 'INS']
    if add and add not in exon_types:
        raise ValueError('Invalid exon type to add, expected %s' % exon_types)
    
    if not all([p in params.keys() for p in ['n_exons', 'ins_range', 'out_prefix', 'block_range']]):
        raise ValueError('Some parameters missing')
    
    n_exons = params['n_exons']
    block_range = params['block_range']
    ins_range = params['ins_range']
    control_fasta = '%s-control.fasta' % params['out_prefix']
    case_fasta = '%s-case.fasta' % params['out_prefix']
    
    # get sequence for tx1    
    seq1, strand1, ex1_list = get_tx_seq(tx1, all_exons, genome_fasta, n_exons, control_fasta)
    
    # add to fusion list
    fusion_parts = [tx1]

    # extended or novel exon
    ext_seq, bloc = '', ''
    if add == 'EE':
        ext_seq, bloc = get_exon_seq(ex1_list, strand1, all_exons, 
                                     genome_fasta, block_range)
    elif add == 'NE':
        ext_seq, bloc = get_exon_seq(ex1_list, strand1, all_exons, 
                                     genome_fasta, block_range, extended=False)
    elif add == 'INS':
        ext_seq = get_random_seq(ins_range)
        bloc = ext_seq
    fusion_parts.append(bloc)

    seq2 = ''
    if tx2:
        # get sequence for tx1
        seq2, strand2, ex2_list = get_tx_seq(tx2, all_exons, genome_fasta, n_exons,
                                             control_fasta, front=False)

        # add to fusion list
        fusion_parts.append(tx2)
    else:
        # unpartnered fusion
        block_seq, seed = get_random_block(all_exons, gene_trees, genome_fasta, block_range)
        seq2 = [s for s in block_seq.values()]
        bloc = ''.join([k for k in block_seq.keys()])
        fusion_parts.append(bloc)

    seq = ''.join(seq1 + [ext_seq] + seq2)
    name = '%s:%s:%s' % (tx1, bloc, tx2) if tx2 else '%s:%s' % (tx1, bloc)
    
    # write output
    with open(case_fasta, 'a') as fout:
        fout.write('>%s\n' % name)
        fout.write(seq + '\n')
        
    return fusion_parts