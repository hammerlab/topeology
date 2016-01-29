# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pandas as pd
import numpy as np
from pepdata import pmbec
from statsmodels.stats.moment_helpers import cov2corr
from skbio.alignment import StripedSmithWaterman
from collections import defaultdict
import itertools
import contextlib
import sys
from six import StringIO

from .iedb_data import get_iedb_epitopes

AMINO_ACID_LETTERS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                      'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                      'T', 'W', 'Y', 'V']
INVALID_AMINO_ACID_LETTERS = ['B', 'Z', 'X', '*']

# Taken from http://stackoverflow.com/questions/2828953
@contextlib.contextmanager
def no_stdout():
    save_stdout = sys.stdout
    sys.stdout = StringIO()
    yield
    sys.stdout = save_stdout

def get_pmbec():
    """Return a PMBEC correlation 2D dictionary."""
    # Silence stdout, since read_coefficients prints to stdout
    # TODO: Just fix pepdata.pmbec to not do this.
    with no_stdout():
        pmbec_coeffs = pmbec.read_coefficients()
        pmbec_coeffs_df = pd.DataFrame(pmbec_coeffs)

    # Use correlation rather than covariance
    pmbec_df = pd.DataFrame(cov2corr(pmbec_coeffs_df))
    pmbec_df.index = pmbec_coeffs_df.index
    pmbec_df.columns = pmbec_coeffs_df.columns

    # Include invalid letters, as Smith-Waterman expects substitution matrix values for them
    pmbec_dict = defaultdict(dict)
    pmbec_dict.update(pmbec_df.to_dict())
    valid_letters = set(pmbec_dict.keys())
    all_letters = valid_letters.union(INVALID_AMINO_ACID_LETTERS)
    for letter_i in all_letters:
        for letter_j in all_letters:
            if not(letter_i in valid_letters and letter_j in valid_letters):
                # We dont need lower than 0, as Smith-Waterman sets negative scores to 0
                pmbec_dict[letter_i][letter_j] = 0

    return pmbec_dict

def pmbec_matrix():
    column_order = []
    column_order.extend(AMINO_ACID_LETTERS)
    column_order.extend(INVALID_AMINO_ACID_LETTERS)

    pmbec_df = pd.DataFrame(get_pmbec())

    # Order the matrix according to this column order in both directions
    # (horizontal and vertical amino acids)
    pmbec_df = pmbec_df[column_order]
    pmbec_df = pmbec_df.T
    pmbec_df = pmbec_df[column_order]
    pmbec_df = pmbec_df.T

    return [int(round(val * 100)) for val in pmbec_df.as_matrix().flatten()]

def pmbec_min():
    return pmbec_matrix().min()

def pmbec_max():
    return pmbec_matrix().max()

def matrix_values_apply_func(matrix_dict, func):
    """Apply func to the values of a 2D dictionary."""
    matrix_values = [inner.values() for inner in matrix_dict.values()]
    return func(itertools.chain(*matrix_values))

def similarity_score(seq_a, seq_b, substitution_dict, gap_penalty):
    """
    Use Smith-Waterman to align seq_a and seq_b using the given substitution
    dictionary (2D) and gap penalty (which is used for both gap opening and gap
    extension.
    """
    # StripedSmithWaterman expects str vs. unicode
    seq_a = str(trim_seq(seq_a))
    seq_b = str(trim_seq(seq_b))

    query = StripedSmithWaterman(
        seq_a, protein=True,
        gap_open_penalty=gap_penalty, gap_extend_penalty=gap_penalty,
        substitution_matrix=dict(substitution_dict))
    return query(seq_b)['optimal_alignment_score']

def trim_seq(seq):
    return seq[2:-1]

def get_neoepitopes(epitope_file_path, epitope_lengths):
    """
    Expected header format: sample, epitope
    """
    df_neoepitopes = pd.read_csv(epitope_file_path, dtype=object, header=0)

    # Only certain lengths
    df_neoepitopes['epitope_length'] = df_neoepitopes['epitope'].apply(len)
    df_neoepitopes = df_neoepitopes[df_neoepitopes['epitope_length'].isin(epitope_lengths)]

    df_neoepitopes.reset_index(drop=True, inplace=True)
    return df_neoepitopes

def multiply_and_round_dict(dict_2d, scalar):
    """
    Multiply values in a 2D dictionary by a scalar, and round the resultant values
    to the nearest integer.
    """
    new_dict = defaultdict(dict)
    for key_i in dict_2d.keys():
        for key_j in dict_2d[key_i].keys():
            new_dict[key_i][key_j] = round(dict_2d[key_i][key_j] * scalar)
    return new_dict

def get_joined_epitopes(epitope_file_path, epitope_lengths):
    df_neoepitopes = get_neoepitopes(epitope_file_path=epitope_file_path,
                                     epitope_lengths=epitope_lengths)
    df_iedb_epitopes = get_iedb_epitopes(epitope_lengths=epitope_lengths)
    return df_neoepitopes.merge(df_iedb_epitopes, on='epitope_length')

def calculate_similarity_from_df(df):
    """
    Given a DataFrame with epitope and iedb_epitope columns, calculate
    a score for every row.
    """
    import imp
    faster = False
    try:
        imp.find_module('pmbecalign')
        faster = True
    except ImportError:
        pass

    if faster:
        from pmbecalign import pmbec_init, pmbec_score
        pmbec_init(pmbec_matrix())
        df['score'] = df.apply(
            lambda row: pmbec_score(trim_seq(row['epitope']),
                                    trim_seq(row['iedb_epitope'])),
            axis=1)
    else:
        # Multiply by 100 to get integers, as StripedSmithWaterman expects integers
        pmbec_dict = get_pmbec()
        multiply_scalar = 100.0
        pmbec_dict = multiply_and_round_dict(pmbec_dict, multiply_scalar)
        pmbec_min = abs(matrix_values_apply_func(pmbec_dict, min))
        df['score'] = df.apply(
            lambda row: similarity_score(row['epitope'], row['iedb_epitope'],
                                         substitution_dict=pmbec_dict,
                                         gap_penalty=pmbec_min), axis=1)
        df.score = df.score.apply(lambda score: float(score) / multiply_scalar)

    return df

def compare(epitope_file_path, epitope_lengths=[8, 9, 10, 11]):
    """
    Given a neoepitope file path, compare each epitope with IEDB and
    return a DataFrame with resultant scores.

    Output columns: sample, epitope, iedb_epitope, score
    """
    df_joined = get_joined_epitopes(epitope_file_path=epitope_file_path,
                                    epitope_lengths=epitope_lengths)
    return calculate_similarity_from_df(df_joined)[[
        'sample', 'epitope', 'iedb_epitope', 'score']]
