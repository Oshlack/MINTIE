import pytest
import pandas as pd
import refine_annotations as ra
from intervaltree import Interval, IntervalTree

# dummy args
args = type('argparse.Namespace', (object,),
            {'minClip': 20, 'minGap': 7})()
ra.set_globals(args)

@pytest.mark.parametrize('coord,expected', [((150, 160), True),
                                            ((90, 101), False),
                                            ((198, 205), False),
                                            ((193, 200), True),
                                            ((100, 107), True)])
def test_check_overlap_del(coord, expected):
    ex_trees = {}
    ref_tree = IntervalTree()
    coords = [(100, 200)]
    for s,e in coords:
        ref_tree.addi(s, e)
    ex_trees['chr1'] = ref_tree

    s, e = coord
    assert ra.check_overlap(ex_trees, 'chr1', s, e, True) == expected

@pytest.mark.parametrize('coord,expected', [((150, 160), True),
                                            ((90, 150), False),
                                            ((150, 250), False),
                                            ((301, 399), True),
                                            ((150, 350), False)])
def test_overlaps_same_exon(coord, expected):
    ex_trees = {}
    ref_tree = IntervalTree()
    coords = [(100, 200),
              (300, 400)]
    for s,e in coords:
        ref_tree.addi(s, e)
    ex_trees['chr1'] = ref_tree

    s, e = coord
    pos1 = 'chr1:%s(+)' % s
    pos2 = 'chr1:%s(-)' % e
    sv = pd.Series({'pos1': pos1, 'pos2': pos2})
    assert ra.overlaps_same_exon(sv, ex_trees) == expected
