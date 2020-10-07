import pytest
import pandas as pd
import numpy as np
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
coords = [(100, 201),
          (300, 401)]
for s,e in coords:
    ref_tree.addi(s, e)
ex_trees['chr1'] = ref_tree

# make reference -- chrom2
ref_tree = IntervalTree()
coords = [(500, 601),
          (800, 901)]
for s,e in coords:
    ref_tree.addi(s, e)
ex_trees['chr2']= ref_tree

@pytest.mark.parametrize('params,expected', [(['AG', 0, '+'], 0),
                                             (['AC', 0, '-'], 0),
                                             (['AC', 0, '+'], 1),
                                             (['AC', 1, '-'], 2),
                                             (['AC', 1, '+'], 2),
                                             (['GT', 1, '+'], 0),
                                             (['GT', 0, '+'], 2),
                                             (['GT', 1, '-'], 1),
                                             (['TG', 1, '+'], 2),
                                             (['TG', 0, '+'], 1),
                                             (['TG', 1, '-'], 2),
                                             (['TG', 0, '-'], 2)])
def test_get_diff_count(params, expected):
    motif, side, strand = params
    sense = strand == '+'
    diff  = ra.get_diff_count(motif, side = side, sense = sense)
    assert diff == expected

@pytest.mark.parametrize('params,expected', [(['AG,GT', 0], True),
                                            (['AC,CT', 0], True),
                                            (['CA,GT', 0], False),
                                            (['AG,GT', 1], True),
                                            (['AG,CT', 1], True),
                                            (['AG,CC', 1], False),
                                            (['AG,CC', 2], True),
                                            (['TT,CC', 2], False),
                                            (['AG,CC', 3], True),
                                            ([',', 0], False),
                                            ([',', 2], False),
                                            (['AG,', 0], True),
                                            (['AC,', 0], True),
                                            ([',GT', 0], True),
                                            ([',CT', 0], True),
                                            (['CA,', 0], False),
                                            (['CA,', 1], False),
                                            (['CA,', 2], False),
                                            (['CA,', 3], True),
                                            ([',CA', 1], True),
                                            ([',CA', 2], True),
                                            ([',TC', 2], False),
                                            ([',AC', 0], False)])
def test_check_valid_motif(params, expected):
    seq, mismatches = params
    block_seqs = seq.split(',')
    left_idx = '' if block_seqs[0] == '' else 0
    right_idx = '' if block_seqs[1] == '' else 1
    valid, motif = ra.check_valid_motif(left_idx, right_idx, block_seqs, mismatches)
    assert valid == expected and motif == ''.join(block_seqs)

@pytest.mark.parametrize('coord,expected', [((150, 160), 10),
                                            ((90, 100), 0),
                                            ((90, 101), 1),
                                            ((198, 205), 3),
                                            ((193, 200), 7),
                                            ((201, 208), 0),
                                            ((100, 107), 7),
                                            ((50, 90), float('nan'))])
def test_get_overlap_size(coord, expected):
    s, e = coord
    if np.isnan(expected):
        assert np.isnan(ra.get_overlap_size(ex_trees, 'chr1', s, e))
    else:
        assert ra.get_overlap_size(ex_trees, 'chr1', s, e) == expected

@pytest.mark.parametrize('coord,expected', [((150, 160), True),
                                            ((90, 101), False),
                                            ((198, 205), False),
                                            ((193, 200), True),
                                            ((100, 107), True)])
def test_check_overlap_with_size(coord, expected):
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


@pytest.mark.parametrize('coord,expected', [((210, 250, 'NE'), False),
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
