import pytest
import pandas as pd
import tempfile
import create_ec_count_matrix as ecm
import numpy as np
from collections import OrderedDict

s1_eq_classes = '3\n' \
                 + '4\n' \
                 + 'TX1\n' \
                 + 'TX2\n' \
                 + 'TX3\n' \
                 + '3\t0\t1\t2\t16\n' \
                 + '2\t0\t1\t6\n' \
                 + '2\t1\t2\t2\n' \
                 + '1\t0\t6\n'
s2_eq_classes = '3\n' \
                 + '4\n' \
                 + 'TX1\n' \
                 + 'TX2\n' \
                 + 'TX3\n' \
                 + '3\t0\t1\t2\t19\n' \
                 + '2\t0\t1\t2\n' \
                 + '2\t1\t2\t8\n' \
                 + '1\t0\t2\n'

output_matrix = 's1\ts2\n' \
                + '0\t6\t2\n' \
                + '0|1\t6\t2\n' \
                + '0|1|2\t16\t19\n' \
                + '1|2\t2\t8'

output_df = 's1\ts2\tec_names\ttx_ids\ttranscript\n' \
                + '0\t6\t2\tec1\t0\tTX1\n' \
                + '0|1\t6\t2\tec2\t0\tTX1\n' \
                + '0|1\t6\t2\tec2\t1\tTX2\n' \
                + '0|1|2\t16\t19\tec3\t0\tTX1\n' \
                + '0|1|2\t16\t19\tec3\t1\tTX2\n' \
                + '0|1|2\t16\t19\tec3\t2\tTX3\n' \
                + '1|2\t2\t8\tec4\t1\tTX2\n' \
                + '1|2\t2\t8\tec4\t2\tTX3\n'

def test_load_ecs():
    with tempfile.NamedTemporaryFile() as ec_file:
        ec_file.write(bytes(s1_eq_classes, 'utf-8'))
        ec_file.flush()

        ec_txs, ec_counts = ecm.load_ecs(ec_file.name)
        assert ec_txs == ['0|1|2', '0|1', '1|2', '0']
        assert ec_counts == [16, 6, 2, 6]

def test_build_ec_matrix():
    with tempfile.NamedTemporaryFile() as ec_file:
        ec_file.write(bytes(s1_eq_classes, 'utf-8'))
        ec_file.flush()
        tx_lookup = ecm.get_tx_lookup(ec_file.name)
        ec1 = ecm.load_ecs(ec_file.name)
    
    with tempfile.NamedTemporaryFile() as ec_file:
        ec_file.write(bytes(s2_eq_classes, 'utf-8'))
        ec_file.flush()
        ec2 = ecm.load_ecs(ec_file.name)

    sample_ecs = [ec1, ec2]
    sample_names = ['s1', 's2']
    with tempfile.NamedTemporaryFile() as matrix_file:
        matrix_file.write(bytes(output_matrix, 'utf-8'))
        matrix_file.flush()
        out_matrix = pd.read_csv(matrix_file.name, sep='\t')

    ec_matrix = ecm.build_ec_matrix(sample_ecs, sample_names)
    assert np.all(out_matrix == ec_matrix)
    
    with tempfile.NamedTemporaryFile() as df_file:
        df_file.write(bytes(output_df, 'utf-8'))
        df_file.flush()
        out_df = pd.read_csv(df_file.name, sep='\t')
        out_df['tx_ids'] = out_df.tx_ids.map(str)
    
    ec_matrix = ecm.construct_dataframe(ec_matrix, tx_lookup)
    assert np.all(out_df == ec_matrix)

