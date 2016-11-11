import nose
from nose.tools import eq_, ok_
import pandas as pd

from topeology.calculate import compare, calculate_similarity_from_df
from .data import data_path

TEST_EPITOPES = data_path('test_epitopes.csv')

def test_single_comparison():
    # TODO: Figure out why this doesn't work when trim_seq isn't used
    # Error: https://github.com/ekg/vcflib/blob/master/src/ssw.c#L556
    df = pd.DataFrame({'Sample': '001',
                       'epitope': ['AAALPGKCGV'],
                       'iedb_epitope': ['EFKEFAAGRR']})
    df_scores = calculate_similarity_from_df(df)

    # TODO: This score is currently unverified
    eq_(df_scores.score.max(), 2.38)

def test_single_comparison_ignore_seqalign():
    df = pd.DataFrame({'Sample': '001',
                       'epitope': ['AAALPGKCGV'],
                       'iedb_epitope': ['EFKEFAAGRR']})
    df_scores = calculate_similarity_from_df(df, ignore_seqalign=True)

    # TODO: This score is currently unverified
    eq_(df_scores.score.max(), 2.38)

def test_calculate_similarity():
    df_scores = compare(epitope_file_path=TEST_EPITOPES,
                        epitope_lengths=[8, 9, 10, 11])
    ok_(df_scores.score.mean() > 1.0)
