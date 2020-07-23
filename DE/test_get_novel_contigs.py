import pytest
import tempfile
import pandas as pd
import numpy as np
import get_novel_contigs as gnc
from collections import OrderedDict

tx_ref = ['>tx1 GENE\n', 'AGTA\n']
ecm_dict = OrderedDict([('sample', [5, 5, 10]),
                        ('ec_names', ['ec1', 'ec1', 'ec2']),
                        ('transcript', ['tx1', 'k79_1', 'k49_1'])])
ecm = pd.DataFrame.from_dict(ecm_dict)

def test_parse_args():
    args = gnc.parse_args(['ec_matrix_file.txt',
                           'tx_ref.fasta',
                           'sample_denovo_filt.fasta'])
    assert args.ecm_file == 'ec_matrix_file.txt'
    assert args.ref_tx_fasta == 'tx_ref.fasta'
    assert args.denovo_fasta == 'sample_denovo_filt.fasta'

def test_get_ref_txs():
    with tempfile.NamedTemporaryFile() as fasta:
        for item in tx_ref:
            fasta.write(bytes(item, 'utf-8'))
            fasta.flush()
        assert gnc.get_ref_txs(fasta.name) == ['tx1']

def test_get_tx_ec():
    res_dict = OrderedDict([('ec_names', ['ec2']),
                            ('contig', ['k49_1']),
                            ('n_contigs_in_ec', [1]),
                            ('contigs_in_EC', ['k49_1']),
                            ('case_reads', [10])])
    results = pd.DataFrame.from_dict(res_dict)
    assert np.all(gnc.get_tx_ec(ecm, ['tx1'], '/tmp').reset_index(drop = True) == results)
