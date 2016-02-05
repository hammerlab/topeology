import nose
from nose.tools import eq_, ok_
import pandas as pd

from topeology.calculate import compare, calculate_similarity_from_df
from .data import data_path

TEST_EPITOPES = data_path('test_epitopes.csv')

def test_fast_versus_slow():
    df = pd.DataFrame({'Sample': '001',
                       'epitope': ['AAAAL'] * 1000,
                       'iedb_epitope': ['AGGGT'] * 1000})
    from timeit import Timer
    slow_timer = Timer(lambda: calculate_similarity_from_df(df, ignore_seqalign=True))
    slow_time = slow_timer.timeit(number=3)
    fast_timer = Timer(lambda: calculate_similarity_from_df(df, ignore_seqalign=False))
    fast_time = fast_timer.timeit(number=3)
    ok_(slow_time > fast_time * 5,
        ("The fast similarity calculation (%0.2f) should be at least "
         "5x faster than the slower version (%0.2f).") % (
             fast_time, slow_time))
