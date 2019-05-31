import pytest
import pandas as pd
import annotate_contigs as ac
import refine_annotations as ra
from intervaltree import Interval, IntervalTree

@pytest.mark.parametrize('letter,next_letter', [('a', 'b'),
                                                ('g', 'h'),
                                                ('z', 'A')])
def test_get_next_letter(letter, next_letter):
    assert ac.get_next_letter(letter) == next_letter

@pytest.mark.parametrize('letter,exception', [('Z', IndexError),
                                              ('1', IndexError)])
def test_get_next_letter_exception(letter, exception):
    with pytest.raises(exception):
        ac.get_next_letter(letter)

def test_get_next_id():
    assert ac.get_next_id('k49_123') == 'k49_123a'
    assert ac.get_next_id('k49_123') == 'k49_123b'

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
