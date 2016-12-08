import nose
from nose.tools import eq_, ok_
import pandas as pd

from topeology.calculate import compare, calculate_similarity_from_df
from .data import data_path

TEST_EPITOPES = data_path('test_epitopes.csv')
TEST_EPITOPES_WILDTYPE = data_path('test_epitopes_wildtype.csv')

def test_single_comparison():
    # TODO: Figure out why this doesn't work when trim_seq isn't used
    # Error: https://github.com/ekg/vcflib/blob/master/src/ssw.c#L556
    df = pd.DataFrame({'sample': '001',
                       'epitope': ['AAALPGKCGV'],
                       'iedb_epitope': ['EFKEFAAGRR']})
    df_scores = calculate_similarity_from_df(df)

    # TODO: This score is currently unverified
    eq_(df_scores.score.max(), 2.38)

def test_single_comparison_ignore_seqalign():
    df = pd.DataFrame({'sample': '001',
                       'epitope': ['AAALPGKCGV'],
                       'iedb_epitope': ['EFKEFAAGRR']})
    df_scores = calculate_similarity_from_df(df, ignore_seqalign=True)

    # TODO: This score is currently unverified
    eq_(df_scores.score.max(), 2.38)

def test_calculate_similarity():
    df_scores = compare(epitopes=TEST_EPITOPES,
                        epitope_lengths=[8, 9, 10, 11])
    ok_(df_scores.score.mean() > 1.0)

def test_wildtype():
    df_scores = compare(epitopes=TEST_EPITOPES_WILDTYPE,
                        epitope_lengths=[8, 9, 10, 11],
                        include_wildtype=True)
    ok_(df_scores.score.mean() > 1.0)
    ok_(df_scores.score_wt.mean() > 1.0)

def test_no_negative_scores():
    """
    AAFFFLVVL has resulted in negative scores in the past.
    """
    df = pd.DataFrame({'sample': '001',
                       'epitope': ['AAFFFLVVL'] * 8,
                       'iedb_epitope': ['DPRRRSRNL',
                                        'GLDRNSGNY',
                                        'IMNRRKRSV',
                                        'KPNRNGGGY',
                                        'KTDNNNSNF',
                                        'LLHDRQHSI',
                                        'SLKKNSRSL',
                                        'EFKEFAAGRR']})
    df_scores_no_seqalign = calculate_similarity_from_df(df, ignore_seqalign=True)
    df_scores = calculate_similarity_from_df(df, ignore_seqalign=False)
    ok_(df_scores.score.min() >= 0)
    ok_(df_scores.equals(df_scores_no_seqalign))
