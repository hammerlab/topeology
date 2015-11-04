import nose
from nose.tools import eq_, ok_
from topeology import calculate_similarity
from .data import data_path

TEST_EPITOPES = data_path('test_epitopes.csv')

def test_calculate_similarity():
    df_scores = calculate_similarity(TEST_EPITOPES)
    assert df_scores.score.mean() == 1
