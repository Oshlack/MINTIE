import pytest
import pandas as pd
import numpy as np
import make_supertranscript as ms


@pytest.mark.parametrize('seq,expected', [('AGTC', 'GACT'),
                                          ('NAGTC', 'GACTN'),
                                          ('A', 'T'),
                                          ('C', 'G'),
                                          ('T', 'A'),
                                          ('G', 'C')])
def test_reverse_complement(seq, expected):
    assert ms.reverse_complement(seq) == expected

def test_get_contig_genes():
    con_info = {'overlapping_genes': ['A:B', 'B:', 'C']}
    con_info = pd.DataFrame.from_dict(con_info)
    assert list(ms.get_contig_genes(con_info)) == ['A', 'B']
