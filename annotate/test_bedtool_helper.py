import pytest
import bedtool_helper as bh
import numpy as np
import pandas as pd
import tempfile
from pybedtools import BedTool

# test data
gtf_data = {'chrom': ['1', '1', '1'],
       'source': ['test', 'test', 'test'],
       'feature': ['gene', 'exon', 'exon'],
       'start': [100, 100, 200],
       'end': [250, 150, 250],
       'score': ['.', '.', '.'],
       'strand': ['+', '+', '+'],
       'frame': ['.', '.', '.'],
       'attribute': ['', '', '']}

def test_subset_featuretypes():
    gtf = pd.DataFrame.from_dict(gtf_data)
    g = BedTool.from_dataframe(gtf)
    exons = BedTool(bh.subset_featuretypes(g, 'exon'))
    assert [e[2] for e in exons] == ['exon', 'exon']

def test_add_strand():
    bed = {'chrom': ['1', '1'],
           'start': [100, 200],
           'end': [150, 250]}
    bed = pd.DataFrame.from_dict(bed)
    b = BedTool.from_dataframe(bed).remove_invalid().sort().merge()
    bex = b.each(bh.add_strand, '+')
    assert [x.strand for x in bex] == ['+', '+']

def test_get_block_seqs():
    with tempfile.NamedTemporaryFile() as fa_tmp:
        seq = np.random.choice(list('AGTC'), 400)
        seq = ''.join(seq)
        fa_tmp.write(bytes('>1\n', 'utf-8'))
        fa_tmp.write(bytes(seq, 'utf-8'))
        fa_tmp.flush()

        gtf = pd.DataFrame.from_dict(gtf_data)
        g = BedTool.from_dataframe(gtf)
        g = g.sequence(fi=fa_tmp.name, s=True)
        block_seqs = bh.get_block_seqs(g)
        assert len(block_seqs) == 3
