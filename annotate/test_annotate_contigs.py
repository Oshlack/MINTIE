import pytest
import pandas as pd
import annotate_contigs as ac

@pytest.mark.parametrize('letter,next_letter', [('a', 'b'),
                                                ('g', 'h'),
                                                ('z', 'A'),
                                                ('Z', 'Za'),
                                                ('ZZ', 'ZZa'),
                                                ('Za', 'Zb')])
def test_get_next_letter(letter, next_letter):
    assert ac.get_next_letter(letter) == next_letter

def test_get_next_id():
    assert ac.get_next_id('k49_123') == 'k49_123a'
    assert ac.get_next_id('k49_123') == 'k49_123b'

