import pytest
import bedtool_helper as bh
import pandas as pd
from pybedtools import BedTool

def test_subset_featuretypes():
    gtf = {'chrom': ['1', '1', '1'],
           'source': ['test', 'test', 'test'],
           'feature': ['gene', 'exon', 'exon'],
           'start': [100, 100, 200],
           'end': [250, 150, 250],
           'score': ['.', '.', '.'],
           'strand': ['+', '+', '+'],
           'frame': ['.', '.', '.'],
           'attribute': ['', '', '']}
    gtf = pd.DataFrame.from_dict(gtf)
    g = BedTool.from_dataframe(gtf)
    exons = BedTool(bh.subset_featuretypes(g, 'exon'))
    assert [e[2] for e in exons] == ['exon', 'exon']
