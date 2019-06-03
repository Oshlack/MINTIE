import pytest
import pandas as pd
import refine_annotations as ra
from intervaltree import Interval, IntervalTree


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
