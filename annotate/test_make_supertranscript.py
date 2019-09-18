import pytest
import make_supertranscript as ms


@pytest.mark.parametrize('seq,expected', [('AGTC', 'GACT'),
                                          ('NAGTC', 'GACTN'),
                                          ('A', 'T'),
                                          ('C', 'G'),
                                          ('T', 'A'),
                                          ('G', 'C')])
def test_reverse_complement(seq, expected):
    assert ms.reverse_complement(seq) == expected
