from pepdata import pmbec
import pandas as pd
from collections import defaultdict
from statsmodels.stats.moment_helpers import cov2corr

SPECIAL_LETTERS = ['B', 'Z', 'X', '*']

def pmbec_df(include_special_letters=True):
    pmbec_coeffs = pmbec.read_coefficients()
    pmbec_coeffs_df = pd.DataFrame(pmbec_coeffs)

    # Use correlation rather than covariance
    pmbec_df = pd.DataFrame(cov2corr(pmbec_coeffs_df))
    pmbec_df.index = pmbec_coeffs_df.index
    pmbec_df.columns = pmbec_coeffs_df.columns

    pmbec_dict = pmbec_df.to_dict()

    special_letters = SPECIAL_LETTERS if include_special_letters else []
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

def pmbec_matrix(include_special_letters=True):
    column_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                    'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                    'T', 'W', 'Y', 'V']
    df = pmbec_df(include_special_letters=include_special_letters)
    if include_special_letters:
        column_order.extend(SPECIAL_LETTERS)
    df = df[column_order]
    df = df.T
    df = df[column_order]
    df = df.T

    return df.as_matrix().flatten()

def pmbec_min():
    return pmbec_matrix(include_special_letters=False).min()
    
def pmbec_max():
    return pmbec_matrix(include_special_letters=False).max()
