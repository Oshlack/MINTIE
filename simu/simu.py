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
import itertools
from IPython.core.debugger import set_trace
from intervaltree import Interval, IntervalTree
from pybedtools import BedTool
from Bio import SeqIO

BASES = list('GCAT')
BASE_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
SEED_INIT = 123
MAX_BP_FROM_BOUNDARY = 10 #to place variant for indels

#=====================================================================================================
# Util funcs
#=====================================================================================================

def increment_seed(seed, amount=1):
    seed += amount
    np.random.seed(seed)
    return seed

def reverse_complement(seq):
    if seq == '':
        return ''
    if type(seq) == float and math.isnan(seq):
        return ''
    seq = seq[::-1]
    seq = ''.join([BASE_COMPLEMENT[base.upper()] for base in list(seq)])
    return(seq)

def get_pos_parts(loc):
    loc = loc.split(':')
    chrom = loc[0]
    pos = loc[1].split('(')[0].split('-')
    start, end = [int(p) for p in pos]
    strand = re.search(r'\(([-+])\)', loc[1]).group(1)
    return chrom, start, end, strand

#=====================================================================================================
# Get funcs
#=====================================================================================================

def pick_genes(n, genes_list):
    '''
    Randomly select n genes from given gene list,
    and remove those genes from the genes list.
    '''
    pick_genes = np.random.choice(genes_list, n, replace=False)
    genes_list = [gene for gene in genes_list if gene not in pick_genes]
    return pick_genes, genes_list

def get_tx_chrom(tx, all_exons):
    chrom = [ex.chrom for ex in all_exons if ex['transcript_id'] == tx][0]
    return chrom

def get_transcripts(gene, all_exons, valid_txs=[]):
    '''
    Pick a transcript given from a specified gene.
    By default always picks the first transcript found.
    Optionally, the transcript is checked against a
    valid_txs list.
    '''
    txs = [ex['transcript_id'] for ex in all_exons if get_gene_name(ex) == gene]
    if len(valid_txs) > 0:
        txs = [tx for tx in txs if tx in valid_txs]
    return np.unique(txs)

def get_gene_name(row):
    '''
    Prevents KeyError if gene name missing
    '''
    try:
        return row.attrs['gene_name']
    except KeyError:
        return ''

def get_gene_loc(chrom, gene_trees, gene):
    '''
    Return coordinates of gene
    '''
    lookup = pd.DataFrame(gene_trees[chrom])
    record = lookup[lookup.data==gene]
    loc = '%s:%s-%s' % (chrom, record.begin.values[0], record.end.values[0])
    return loc

def get_valid_txs(all_exons, min_exons):
    '''
    Returns all valid transcripts (with at least [min_exons]
    exons, for fusions or splice variants),as well as all
    genes associated with these valids transcripts.
    '''
    all_txs = [(tx['transcript_id'], get_gene_name(tx)) for tx in all_exons]
    valid_txs = pd.DataFrame(pd.Series(all_txs).value_counts(), columns=['exon_count'])
    valid_txs = valid_txs[valid_txs.exon_count >= min_exons]
    valid_txs = valid_txs.index.values

    all_genes = np.unique([gene for tx, gene in valid_txs if gene != ''])

    return valid_txs, all_genes

def get_chrom_features(chrom, gr):
    '''
    Get all merged intervals on given chromosome
    '''
    chrom_features = gr.filter(lambda x: x.chrom == chrom).merge()
    chrom_features = [(g.start, g.end) for g in gr]

    chrom_tree = IntervalTree()
    [chrom_tree.addi(s, e) for s, e in chrom_features]

    return chrom_tree

#=====================================================================================================
# Junc/gene coordinate funcs
#=====================================================================================================

def get_juncs(tx):
    '''
    Return list of junctions in form
    [(chr, start, end)] from transcript info file
    note that exon *ends* become junction *starts*
    '''
    starts = tx['exonStarts'].split(',')[1:]
    ends = tx['exonEnds'].split(',')[:-1]
    chroms = [tx['chrom']] * len(starts)
    return(list(zip(chroms, ends, starts)))

def build_junc_ref(junc_file):
    '''
    Take junction reference file, generate and
    return dictionary containing all splice
    junction start and end sites
    '''
    genref = pd.read_csv(junc_file, sep='\t')
    junc_info = genref.apply(lambda tx: get_juncs(tx), axis=1)

    juncs = ['%s:%s-%s' % (c, s, e) for jv in junc_info.values for c, s, e in jv]
    junc_dic = {}
    for junc in juncs:
        junc_dic[junc] = True

    return junc_dic

def get_junction(ex1, ex2):
    ex1 = re.split('[:()]', ex1)[:-1]
    ex2 = re.split('[:()]', ex2)[:-1]
    s_a, e_a = [int(e) for e in ex1[1].split('-')]
    s_b, e_b = [int(e) for e in ex2[1].split('-')]

    start, end = e_b, s_a
    if e_a < s_b:
        start, end = e_a, s_b

    assert start < end
    return start, end

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

#=====================================================================================================
# Sequence funcs
#=====================================================================================================

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

def get_intron_seq(ex_list, strand, gr, genome_fasta):
    '''
    Get intron sequence of downstream exon
    '''
    block = ex_list[-1] if strand == '+' else ex_list[0]
    loc = re.compile('[:\-\(\)]').split(block)
    exon_start, exon_end = int(loc[1]), int(loc[2])

    pos = exon_end if strand == '+' else exon_start
    ex = pd.DataFrame(get_chrom_features(loc[0], gr.merge()))

    if strand == '-':
        start, end = max(ex[ex.end < pos].end), pos
    else:
        start, end = pos, min(ex[ex.begin > pos].begin)

    block_bed = '%s\t%d\t%d\t.\t1\t%s' % (loc[0], start, end, strand)
    block_bt = BedTool(block_bed, from_string=True)
    block_seq, bs = get_seq(block_bt, genome_fasta)
    ext_seq = ''.join([bs for bs in block_seq.values()])
    bloc = ''.join([k for k in block_seq.keys()])

    return ext_seq, bloc

def get_random_seq(ins_range):
    '''
    Generate random insertion sequence
    '''
    ins_size = np.random.randint(ins_range[0], ins_range[1])
    ins = np.random.choice(BASES, ins_size)
    ins = ''.join(ins)
    return ins

def get_exon_seq(ex_list, strand, gr, genome_fasta, block_range, extended=True):
    '''
    Extends given exon, or creates a novel downstream
    exon with random size and returns its sequence.
    If a reference exon exists that already extends
    or overlaps the given exon, it will extend or place
    the exon past the overlapping exon.
    '''
    block_size = np.random.randint(block_range[0], block_range[1])
    gap_size = 0 if extended else np.random.randint(block_range[0], block_range[1])

    block = ex_list[-1] if strand == '+' else ex_list[0]
    loc = re.compile('[:\-\(\)]').split(block)
    exon_start, exon_end = int(loc[1]), int(loc[2])

    start = exon_end + gap_size if strand == '+' else exon_start - gap_size - block_size
    end = start + block_size

    # extend exon if overlaps an existing one
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

def get_tx_seq(tx, all_exons, genome_fasta, n_exons=0, front=True, control_fasta=''):
    '''
    Get fusion sequence of given transcript, returning
    sequence of first N exons for transcript 1 (front=True)
    and N exons for transcript 2 (front=False). By default,
    write out wildtype transcript to control reference file.
    '''

    exons = all_exons.filter(lambda x: x['transcript_id'] == tx).saveas()
    tx_seq, s = get_seq(exons, genome_fasta)
    ex_list = [ex for ex in tx_seq.keys()]

    if n_exons > 0:
        # pick N 5' exons for tx1 (front) and N 3' exons for tx2 (back)
        if front:
            ex_list = ex_list[:n_exons] if s == '+' else ex_list[-n_exons:]
        else:
            ex_list = ex_list[-n_exons:] if s == '+' else ex_list[:n_exons]

    # select sequences and reverse order if antisense
    seq = [tx_seq[ex] for ex in ex_list]
    seq = seq if s == '+' else [s for s in reversed(seq)]
    ex_list = ex_list if s == '+' else [ex for ex in reversed(ex_list)]

    if control_fasta != '':
        write_wildtype_sequence(tx_seq, s, control_fasta, tx)

    return seq, s, ex_list

#=====================================================================================================
# Exon and transcript selection
#=====================================================================================================

def pick_valid_exon_tx(gene, all_exons, valid_txs, genome_fasta, block_range, vartype):
    '''
    Pick a transcript that has a valid exon to use for
    extending/creating a novel exon. Return transcript
    name and the valid exon number.
    '''
    vartypes = ['EE', 'NE', 'RI']
    if vartype not in vartypes:
        raise ValueError('Invalid variant type to add, expected %s' % vartypes)

    txs = get_transcripts(gene, all_exons)
    tx = txs[0] # pick first as default
    seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta)

    select = get_exon_to_extend(ex_list, all_exons, block_range, vartype)
    if select is not None:
        return tx, select

    # no valid exons in first transcript selected, we have to go deeper...
    txs = [x for x in txs if x!=tx]
    for tx in txs:
        seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta)
        select = get_exon_to_extend(ex_list, all_exons, block_range, vartype)
        if select is not None:
            break

    return tx, select

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

def get_exon_for_deletion(seq, indel_range):
    '''
    Select suitable exon to delete sequence from
    '''
    select = np.random.randint(len(seq))
    max_delsize = len(seq[select]) - (MAX_BP_FROM_BOUNDARY * 2)
    varmin, varmax = indel_range[0], min(indel_range[1], max_delsize)

    # ensure exon is big enough to delete from, otherwise pick another exon
    if varmin > varmax:
        other_exons = [i for i in range(len(seq)) if i != select]
        for select in other_exons:
            max_delsize = len(seq[select]) - (MAX_BP_FROM_BOUNDARY * 2)
            # ^ max possible variant size including padding from boundary
            varmax = min(indel_range[1], max_delsize)
            if varmin < varmax:
                break
    assert varmin <= varmax # bad transcript if assert fails

    return select, varmin, varmax

def get_exon_for_insertion(seq):
    '''
    Select suitable exon to insert sequence
    '''
    select = np.random.randint(len(seq))

    if len(seq[select]) < (MAX_BP_FROM_BOUNDARY * 2):
        other_exons = [i for i in range(len(seq)) if i != select]
        for select in other_exons:
            if len(seq[select]) >= (MAX_BP_FROM_BOUNDARY * 2): break
    assert len(seq[select]) >= (MAX_BP_FROM_BOUNDARY * 2) # bad transcript if assert fails

    return select

def get_exon_to_extend(ex_list, all_exons, block_range, vartype):
    '''
    Select suitable exon for extending or creating novel
    downstream exon. Criteria:
    1. Must not be last exon.
    2. There are no other extended exons
       that extend past the *selected* exon.
    '''
    loc = re.split('[:\(\)]', ex_list[0])[:-1]
    ex_ref = get_chrom_features(loc[0], all_exons.merge())

    # make sure last exon is not selected and randomise order
    exons = [ex for ex in range(0, len(ex_list)-1)]
    exons = np.random.choice(exons, len(exons), replace=False)

    # max distance from downstream exon
    maxlen = block_range[1] * 2 if vartype == 'NE' else block_range[1]

    # make sure there are no overlaps
    select, olaps = None, True
    for select in exons:
        chrom, start, end, strand = get_pos_parts(ex_list[select])
        if vartype == 'RI':
            # maxlen becomes length to next exon for retained intron vars
            nselect = select
            nchrom, nstart, nend, nstrand = get_pos_parts(ex_list[select + 1])
            maxlen = nstart - end - 1 if strand == '+' else start - nend - 1
        if strand == '+':
            olaps = ex_ref.overlaps(end + 1, end + maxlen)
        else:
            olaps = ex_ref.overlaps(start - maxlen, start - 1)
        if not olaps:
            break

    select = select if not olaps else None
    return select

def truncate_exon(ex, ss, block_range, ex_lookup, right=True):
    '''
    Truncate an exon on the right or left side by block_range size.
    Ensure that exon is not truncated at an existing boundary.
    '''
    max_trunc = min(block_range[1], len(ss) - MAX_BP_FROM_BOUNDARY)
    trunc = np.random.randint(block_range[0], max_trunc)

    chrom, start, end, strand = get_pos_parts(ex)
    trunc_pos = end - trunc if right else start + trunc

    locs = np.concatenate([ex_lookup.begin.values, ex_lookup.end.values])
    while trunc_pos in locs:
        trunc_pos = trunc_pos - 1 if right else start - 1
    assert trunc > MAX_BP_FROM_BOUNDARY

    ss = ss[:(len(ss)-trunc)] if right else ss[trunc:]
    return ss, trunc

#=====================================================================================================
# Write functions
#=====================================================================================================

def write_output(seq, tx, vartype, case_fasta):
    seq = ''.join(seq)
    name = '%s(%s)' % (tx, vartype)
    with open(case_fasta, 'a') as fout:
        fout.write('>%s\n' % name)
        fout.write(seq + '\n')

def write_wildtype_sequence(seq_dict, strand, output_file, name):
    '''
    Writes sequence dictionary to output file
    '''
    seq = [seq_dict[ex] for ex in seq_dict.keys()]
    seq = seq if strand == '+' else [s for s in reversed(seq)]

    with open(output_file, 'a') as fout:
        fout.write('>%s\n' % name)
        fout.write(''.join(seq) + '\n')

#=====================================================================================================
# Make variant functions
#=====================================================================================================

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

    if not all([p in params.keys() for p in ['n_exons', 'ins_range',
                                             'out_prefix', 'block_range']]):
        raise ValueError('Some parameters missing')

    n_exons = params['n_exons']
    block_range = params['block_range']
    ins_range = params['ins_range']
    control_fasta = '%s-control.fasta' % params['out_prefix']
    case_fasta = '%s-case.fasta' % params['out_prefix']

    # get sequence for tx1
    seq1, strand1, ex1_list = get_tx_seq(tx1, all_exons, genome_fasta,
                                         n_exons=n_exons, control_fasta=control_fasta)

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
        seq2, strand2, ex2_list = get_tx_seq(tx2, all_exons, genome_fasta,
                                             n_exons=n_exons, front=False,
                                             control_fasta=control_fasta)
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

def write_indel(tx, all_exons, genome_fasta, indel_range, out_prefix, vartype='DEL'):
    '''
    Write transcript with deletion, insertion or ITD in random exon
    '''
    vartypes = ['DEL', 'INS', 'ITD']
    if vartype not in vartypes:
        raise ValueError('Invalid variant type to add, expected %s' % vartypes)

    control_fasta = '%s-control.fasta' % out_prefix
    case_fasta = '%s-case.fasta' % out_prefix
    seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta, control_fasta=control_fasta)
    varsize, select = None, None

    if vartype == 'DEL':
        select, varmin, varmax = get_exon_for_deletion(seq, indel_range)

        delseq = seq[select]
        varsize = np.random.randint(varmin, varmax)
        delpos = np.random.randint(MAX_BP_FROM_BOUNDARY,
                                   (len(delseq) - varsize - MAX_BP_FROM_BOUNDARY))

        # perform deletion
        delseq = delseq[:delpos] + delseq[(delpos + varsize):]
        seq[select] = delseq
    elif vartype == 'INS':
        select = get_exon_for_insertion(seq)

        select_seq = seq[select]
        maxpos = len(select_seq) - MAX_BP_FROM_BOUNDARY
        ins_pos = np.random.randint(MAX_BP_FROM_BOUNDARY, maxpos)

        # perform insertion
        insertion = get_random_seq(indel_range)
        seq[select] = select_seq[:ins_pos] + insertion + select_seq[ins_pos:]

        varsize = len(insertion)
    elif vartype == 'ITD':
        # here we use the same parameters for placing our insertion as with deletions
        select, varmin, varmax = get_exon_for_deletion(seq, indel_range)

        # make ITD sequence
        select_seq, varsize = seq[select], np.random.randint(varmin, varmax)
        maxpos = len(select_seq) - varsize - MAX_BP_FROM_BOUNDARY
        itd_pos = np.random.randint(MAX_BP_FROM_BOUNDARY, maxpos)

        itd_seq = select_seq[itd_pos:(itd_pos + varsize)]
        seq[select] = select_seq[:itd_pos] + itd_seq + select_seq[itd_pos:]

    write_output(seq, tx, vartype, case_fasta)
    return varsize, select+1

def write_large_tsv(tx, all_exons, genome_fasta, out_prefix, exons_range, vartype='PTD'):
    '''
    Write transcript with deletion, insertion or ITD in random exon
    '''
    vartypes = ['PTD', 'INV']
    if vartype not in vartypes:
        raise ValueError('Invalid variant type to add, expected %s' % vartypes)

    control_fasta = '%s-control.fasta' % out_prefix
    case_fasta = '%s-case.fasta' % out_prefix
    seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta, control_fasta=control_fasta)

    min_exons, max_exons = exons_range
    max_exons = min(max_exons, len(seq))
    assert min_exons < max_exons <= len(seq)
    n_exons = np.random.randint(min_exons, max_exons)

    select = np.random.randint(len(seq) - n_exons)
    var_seq = seq[select:(select + n_exons)]

    if vartype == 'PTD':
        seq = seq[:select] + var_seq + var_seq + seq[select+n_exons:]
    elif vartype == 'INV':
        var_seq = ''.join(var_seq)
        var_seq = reverse_complement(var_seq)
        seq = seq[:select] + [var_seq] + seq[select+n_exons:]

    write_output(seq, tx, vartype, case_fasta)
    return n_exons, select+1

def write_novel_exon(gene, valid_txs, all_exons, genome_fasta,
                     out_prefix, block_range, gene_trees, vartype='EE'):
    '''
    Write extended exon, novel exon and retained intron variants
    '''
    vartypes = ['EE', 'NE', 'RI']
    if vartype not in vartypes:
        raise ValueError('Invalid variant type to add, expected %s' % vartypes)

    control_fasta = '%s-control.fasta' % out_prefix
    case_fasta = '%s-case.fasta' % out_prefix
    maxgap = block_range[1] if vartype == 'NE' else 0

    seq, strand, ex_list = None, None, None
    tx, select = pick_valid_exon_tx(gene, all_exons, valid_txs,
                                    genome_fasta, block_range, vartype)

    if select is None:
        return '', '', {}

    # now that we know we have a valid tx we can write the wildtype out
    seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta, control_fasta=control_fasta)
    select_ex = ex_list[:(select+1)] if strand == '+' else ex_list[select:]
    # ^ need to cut the exon list to get terminal exon to modify for get_exon_seq

    ext_seq = ''
    if vartype in ['EE', 'NE']:
        extended = vartype == 'EE'
        ext_seq, bloc = get_exon_seq(select_ex, strand, all_exons,
                                     genome_fasta, block_range, extended=extended)
    elif vartype == 'RI':
        ext_seq, bloc = get_intron_seq(select_ex, strand, all_exons, genome_fasta)

    seq = seq[:(select+1)] + [ext_seq] + seq[(select+1):]

    write_output(seq, tx, vartype, case_fasta)
    chrom = get_tx_chrom(tx, all_exons)
    loc = get_gene_loc(chrom, gene_trees, gene)
    stats = {'varsize': len(ext_seq), 'exon': select+1}

    return tx, loc, stats

def write_trunc_exons(tx, all_exons, genome_fasta, out_prefix, block_range):
    '''
    Truncate two adjacent exons to create a novel exon junction
    '''
    control_fasta = '%s-control.fasta' % out_prefix
    case_fasta = '%s-case.fasta' % out_prefix
    min_exon_len = block_range[0] + MAX_BP_FROM_BOUNDARY + 1

    seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta, control_fasta=control_fasta)

    s_min = 0 if strand == '+' else 1
    s_max = len(seq)-1 if strand == '+' else len(seq)
    select = np.random.randint(s_min, s_max)
    next_select = select + 1 if strand == '+' else select - 1

    if len(seq[select]) < min_exon_len or len(seq[next_select]) < min_exon_len:
        other_exons = [i for i in range(s_min, s_max) if i != select]
        for select in other_exons:
            next_select = select + 1 if strand == '+' else select - 1
            if len(seq[select]) >= min_exon_len and len(seq[next_select]) >= min_exon_len:
                break
    assert len(seq[select]) >= min_exon_len and len(seq[next_select]) >= min_exon_len

    chrom = [x.chrom for x in all_exons if x['transcript_id'] == tx][0]
    ex_lookup = pd.DataFrame(get_chrom_features(chrom, all_exons.merge()))

    r1, r2 = strand == '+', strand == '-'
    seq[select], trunc1 = truncate_exon(ex_list[select], seq[select], block_range,
                                        ex_lookup, right=r1)
    seq[next_select], trunc2 = truncate_exon(ex_list[next_select], seq[next_select],
                                             block_range, ex_lookup, right=r2)
    write_output(seq, tx, 'NEJ', case_fasta)

    selected_exons = '%d, %d' % (select+1, next_select+1)
    trunc_lens = '%d, %d' % (trunc1, trunc2)
    return trunc_lens, selected_exons

def write_unannot_splice(gene, all_exons, valid_txs, genome_fasta, out_prefix, junc_ref, gene_trees):
    '''
    Connect exons in random order, checking if they correspond to an
    existing splice junction. If not, return a transcript with this
    junction, otherwise return -1 and '' for both return values if
    not unannotated junction could be created.
    '''
    control_fasta = '%s-control.fasta' % out_prefix
    case_fasta = '%s-case.fasta' % out_prefix

    txs = get_transcripts(gene, all_exons, valid_txs=valid_txs)
    select1, select2, size, strand, chrom = None, None, None, None, None
    for tx in txs:
        seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta)
        chrom = get_tx_chrom(tx, all_exons)

        exons = range(len(seq))
        possible_juncs = [ex for ex in itertools.combinations(exons, 2)]

        for select1, select2 in possible_juncs:
            start, end = get_junction(ex_list[select1], ex_list[select2])
            junc = '%s:%d-%d' % (chrom, start, end)
            if junc not in junc_ref:
                size = end - start
                break
        if junc not in junc_ref:
            break

    if junc in junc_ref:
        # all possible splice variants must already exist....
        return '', '', {}

    # now that we know we have a valid tx we can write the wildtype out
    seq, strand, ex_list = get_tx_seq(tx, all_exons, genome_fasta, control_fasta=control_fasta)

    if select2 < select1:
        select1, select2 = select2, select1
    seq = seq[:(select1+1)] + seq[select2:]

    write_output(seq, tx, 'US', case_fasta)
    loc = get_gene_loc(chrom, gene_trees, gene)
    stats = {'varsize': size, 'exon': '%d-%d' % (select1+1, select2+1)}

    return tx, loc, stats
