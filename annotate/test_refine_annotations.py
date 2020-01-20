import pytest
import pandas as pd
import refine_annotations as ra
from intervaltree import Interval, IntervalTree

###############################################################
# dummy arguments
###############################################################
args = type('argparse.Namespace', (object,),
            {'minClip': 20, 'minGap': 7})()
ra.set_globals(args)

###############################################################
# test data
###############################################################

# make reference -- chrom1
# these data could represent exons or genes
ex_trees = {}
ref_tree = IntervalTree()
coords = [(100, 200),
          (300, 400)]
for s,e in coords:
    ref_tree.addi(s, e)
ex_trees['chr1'] = ref_tree

# make reference -- chrom2
ref_tree = IntervalTree()
coords = [(500, 600),
          (800, 900)]
for s,e in coords:
    ref_tree.addi(s, e)
ex_trees['chr2']= ref_tree

@pytest.mark.parametrize('seq,expected', [('AG,GT', True),
                                          ('AC,CT', True),
                                          ('CA,GT', False),
                                          (',', False),
                                          ('AG,', True),
                                          ('AC,', True),
                                          (',GT', True),
                                          (',CT', True),
                                          ('CA,', False),
                                          (',AC', False)])
def test_is_valid_motif(seq, expected):
    block_seqs = seq.split(',')
    left_idx = '' if block_seqs[0] == '' else 0
    right_idx = '' if block_seqs[1] == '' else 1
    assert ra.is_valid_motif(left_idx, right_idx, block_seqs) == expected

@pytest.mark.parametrize('coord,expected', [((150, 160), True),
                                            ((90, 101), False),
                                            ((198, 205), False),
                                            ((193, 200), True),
                                            ((100, 107), True)])
def test_check_overlap_del(coord, expected):
    s, e = coord
    assert ra.check_overlap(ex_trees, 'chr1', s, e, size = args.minGap) == expected

@pytest.mark.parametrize('coord,expected', [((150, 160), True),
                                            ((90, 101), True),
                                            ((198, 205), True),
                                            ((193, 200), True),
                                            ((100, 107), True),
                                            ((50, 90), False)])
def test_check_overlap(coord, expected):
    s, e = coord
    assert ra.check_overlap(ex_trees, 'chr1', s, e) == expected

def test_get_pos_parts():
    assert ra.get_pos_parts('chr1:100(+)') == ('chr1', 100, '+')
    assert ra.get_pos_parts('chr1:100(-)') == ('chr1', 100, '-')

def test_get_varsize():
    sv = {'pos1': 'chr1:100(+)', 'pos2': 'chr1:200(+)'}
    assert ra.get_varsize(sv) == 100

@pytest.mark.parametrize('coord,expected', [((150, 160), True),
                                            ((90, 150), False),
                                            ((150, 250), False),
                                            ((301, 399), True),
                                            ((150, 350), False),
                                            ((50, 60), False)])
def test_overlaps_same_exon(coord, expected):
    s, e = coord
    pos1 = 'chr1:%s(+)' % s
    pos2 = 'chr1:%s(-)' % e
    sv = pd.Series({'pos1': pos1, 'pos2': pos2})
    assert ra.overlaps_same_exon(sv, ex_trees) == expected


@pytest.mark.parametrize('coord,expected', [(('chr1:100(+)', 'chr1:105', 'DEL'), False),
                                            (('chr1:100(+)', 'chr1:110', 'DEL'), True),
                                            (('chr1:100(+)', 'chr1:119', 'NEJ'), False),
                                            (('chr1:100(+)', 'chr1:121', 'NEJ'), True),
                                            (('chr1:100(+)', 'chr1:101', 'INS'), True),
                                            (('chr1:150(+)', 'chr2:250', 'FUS'), True)])
def test_overlaps_exon(coord, expected):
    pos1, pos2, vartype = coord
    sv = {'pos1': pos1, 'pos2': pos2, 'variant_type': vartype}
    assert ra.overlaps_exon(sv, ex_trees) == expected

@pytest.mark.parametrize('coord,expected', [(('chr1:100(+)', 'chr1:200(+)'), True),
                                            (('chr1:101(+)', 'chr2:400(+)'), True),
                                            (('chr1:450(+)', 'chr2:501(+)'), True),
                                            (('chr1:550(+)', 'chr1:650(+)'), False),
                                            (('chr2:550(+)', 'chr2:650(+)'), True),
                                            (('chr2:150(+)', 'chr2:250(+)'), False)])
def test_overlaps_gene(coord, expected):
    row = {'pos1': coord[0], 'pos2': coord[1]}
    assert ra.overlaps_gene(row, ex_trees) == expected

def test_match_splice_juncs():
    contigs = {'pos1': ['chr1:100(+)', 'chr1:200(+)', 'chr1:400(+)'],
               'pos2': ['chr1:200(+)', 'chr1:300(+)', 'chr1:500(+)'],
               'variant_type': ['NE', 'NEJ', 'INS'],
               'variant_id': ['A', 'B', 'C']}
    contigs = pd.DataFrame.from_dict(contigs)
    assert all(ra.match_splice_juncs(contigs) == pd.Series([True, False, False]))


@pytest.mark.parametrize('coord,expected', [((200, 250, 'NE'), False),
                                            ((100, 150, 'DEL'), True),
                                            ((450, 500, 'FUS'), False)])
def test_vars_overlap_exons(coord, expected):
    s, e, t = coord
    contigs = {'pos1': ['chr1:%d(+)' % s],
               'pos2': ['chr1:%d(+)' % e],
               'variant_type': [t],
               'variant_id': ['A']}
    contigs = pd.DataFrame.from_dict(contigs)
    assert ra.vars_overlap_exon(contigs, ex_trees)[0] == expected

@pytest.mark.parametrize('coord,expected', [((140, 160), True),
                                            ((160, 200), True),
                                            ((140, 145), False)])
def test_get_junc_vars(coord, expected):
    s, e = coord
    contigs = {'pos1': ['chr1:%d(+)' % s],
               'pos2': ['chr1:%d(+)' % e],
               'variant_type': ['NEJ'],
               'variant_id': ['A'],
               'overlaps_exon': [True],
               'large_varsize': e - s > args.minClip}
    contigs = pd.DataFrame.from_dict(contigs)
    assert ('A' in ra.get_junc_vars(contigs, ex_trees, args)) == expected

@pytest.mark.parametrize('coord,expected', [((100, 150, 'NE', 50), False),
                                            ((100, 150, 'DEL', 0), True),
                                            ((100, 100, 'INS', 10), True),
                                            ((100, 100, 'UN', 20), True),
                                            ((100, 100, 'UN', 19), False)])
def test_get_tsv_vars(coord, expected):
    s, e, t, cv = coord
    contigs = {'pos1': ['chr1:%d(+)' % s],
               'pos2': ['chr1:%d(+)' % e],
               'variant_type': [t],
               'variant_id': ['A'],
               'contig_varsize': [cv]}
    contigs = pd.DataFrame.from_dict(contigs)
    contigs['varsize'] = [cv] if t == 'UN' else [e - s]
    contigs['overlaps_exon'] = ra.vars_overlap_exon(contigs, ex_trees)
    assert ('A' in ra.get_tsv_vars(contigs)) == expected


def test_get_fusion_vars():
    contigs = {'pos1': ['chr1:100(+)', 'chr1:150(+)', 'chr1:400(+)'],
               'pos2': ['chr1:150(+)', 'chr2:200(+)', 'chr1:500(+)'],
               'variant_type': ['EE', 'FUS', 'DEL'],
               'variant_id': ['A', 'B', 'C'],
               'contig_id': ['C1', 'C1', 'C2']}
    contigs = pd.DataFrame.from_dict(contigs)
    assert all([fv in ['A', 'B'] for fv in ra.get_fusion_vars(contigs)])
