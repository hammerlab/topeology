from pepdata import pmbec
import pandas as pd
from collections import defaultdict
from statsmodels.stats.moment_helpers import cov2corr

def pmbec_df():
    pmbec_coeffs = pmbec.read_coefficients()
    pmbec_coeffs_df = pd.DataFrame(pmbec_coeffs)

    # Use correlation rather than covariance
    pmbec_df = pd.DataFrame(cov2corr(pmbec_coeffs_df))
    pmbec_df.index = pmbec_coeffs_df.index
    pmbec_df.columns = pmbec_coeffs_df.columns

    pmbec_dict = pmbec_df.to_dict()

    special_letters = ['B', 'Z', 'X', '*']
    normal_letters = pmbec_dict.keys()
    all_letters = special_letters + normal_letters
    new_dict = defaultdict(dict)
    for letter in all_letters:
        for other_letter in all_letters:
            new_dict[letter][other_letter] = -100
            if letter in normal_letters and other_letter in normal_letters:
                new_dict[letter][other_letter] = int(round(
                    pmbec_dict[letter][other_letter] * 100))
    return pd.DataFrame(new_dict)














