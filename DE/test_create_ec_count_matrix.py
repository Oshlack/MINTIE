import pytest
import pandas as pd
import tempfile
import create_ec_count_matrix as ecm
import numpy as np
from collections import OrderedDict

# equivalence class input 'files' for 2 samples
S1_EC_FILE = '3\n' \
                 + '4\n' \
                 + 'TX1\n' \
                 + 'TX2\n' \
                 + 'TX3\n' \
                 + '3\t0\t1\t2\t16\n' \
                 + '2\t0\t1\t6\n' \
                 + '2\t1\t2\t2\n' \
                 + '1\t0\t6\n'
S2_EC_FILE = '3\n' \
                 + '4\n' \
                 + 'TX1\n' \
                 + 'TX2\n' \
                 + 'TX3\n' \
                 + '3\t0\t1\t2\t19\n' \
                 + '2\t0\t1\t2\n' \
                 + '2\t1\t2\t8\n' \
                 + '1\t0\t2\n'

# inter mediate output matrix
OUTPUT_MATRIX = 's1\ts2\n' \
                + '0\t6\t2\n' \
                + '0|1\t6\t2\n' \
                + '0|1|2\t16\t19\n' \
                + '1|2\t2\t8'

# final output
OUTPUT_DF = 's1\ts2\tec_names\ttx_ids\ttranscript\n' \
                + '0\t6\t2\tec1\t0\tTX1\n' \
                + '0|1\t6\t2\tec2\t0\tTX1\n' \
                + '0|1\t6\t2\tec2\t1\tTX2\n' \
                + '0|1|2\t16\t19\tec3\t0\tTX1\n' \
                + '0|1|2\t16\t19\tec3\t1\tTX2\n' \
                + '0|1|2\t16\t19\tec3\t2\tTX3\n' \
                + '1|2\t2\t8\tec4\t1\tTX2\n' \
                + '1|2\t2\t8\tec4\t2\tTX3\n'

def make_ec_matrix(tmp_path):
    ec_file1 = tmp_path / 'ec_file1.txt'
    ec_file1.write_text(S1_EC_FILE)

    ec_file2 = tmp_path / 'ec_file2.txt'
    ec_file2.write_text(S2_EC_FILE)

    ec1 = ecm.load_ecs(ec_file1)
    ec2 = ecm.load_ecs(ec_file2)

    sample_ecs = [ec1, ec2]
    sample_names = ['s1', 's2']

    return ecm.build_ec_matrix(sample_ecs, sample_names)

def test_load_ecs(tmp_path):
    ec_file = tmp_path / 'ec_file1.txt'
    ec_file.write_text(S1_EC_FILE)

    ec_txs, ec_counts = ecm.load_ecs(ec_file)
    assert ec_txs == ['0|1|2', '0|1', '1|2', '0']
    assert ec_counts == [16, 6, 2, 6]

def test_build_ec_matrix(tmp_path):
    ec_matrix = make_ec_matrix(tmp_path)

    matrix_file = tmp_path / 'output_matrix.txt'
    matrix_file.write_text(OUTPUT_MATRIX)
    out_matrix = pd.read_csv(matrix_file, sep='\t')

    assert np.all(out_matrix == ec_matrix)

def test_construct_dataframe(tmp_path):
    ec_file = tmp_path / 'ec_file1.txt'
    ec_file.write_text(S1_EC_FILE)
    tx_lookup = ecm.get_tx_lookup(ec_file)

    ec_matrix = make_ec_matrix(tmp_path)
    ec_df = ecm.construct_dataframe(ec_matrix, tx_lookup)

    df_file = tmp_path / 'output_df.txt'
    df_file.write_text(OUTPUT_DF)
    out_df = pd.read_csv(df_file, sep='\t')
    out_df['tx_ids'] = out_df.tx_ids.map(str)

    assert np.all(out_df == ec_df)
