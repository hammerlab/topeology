# Copyright (c) 2016. Mount Sinai School of Medicine
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

import contextlib
import sys
from six import StringIO
from statsmodels.stats.moment_helpers import cov2corr
from collections import defaultdict
from pepdata import pmbec
import pandas as pd

# Taken from http://stackoverflow.com/questions/2828953
@contextlib.contextmanager
def no_stdout():
    save_stdout = sys.stdout
    sys.stdout = StringIO()
    yield
    sys.stdout = save_stdout

AMINO_ACID_LETTERS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                      'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                      'T', 'W', 'Y', 'V']
INVALID_AMINO_ACID_LETTERS = ['B', 'Z', 'X', '*']

# We multiply PMBEC values by 100 to use them in Smith-Waterman alignment
# TODO: Do this more cleanly.
MULTIPLIER = 100.0

class PMBEC(object):
    """Representation of the PMBEC correlation matrix."""

    def __init__(self):
        self.d = self.create_dict()

    def create_dict(self):
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

    def as_int_dict(self, multiplier=MULTIPLIER):
        new_dict = defaultdict(dict)
        for key_i in self.d.keys():
            for key_j in self.d[key_i].keys():
                new_dict[key_i][key_j] = round(self.d[key_i][key_j] * multiplier)
        return new_dict

    def as_int_list(self, multiplier=MULTIPLIER):
        column_order = []
        column_order.extend(AMINO_ACID_LETTERS)
        column_order.extend(INVALID_AMINO_ACID_LETTERS)

        pmbec_df = pd.DataFrame(self.d)

        # Order the matrix according to this column order in both directions
        # (horizontal and vertical amino acids)
        pmbec_df = pmbec_df[column_order]
        pmbec_df = pmbec_df.T
        pmbec_df = pmbec_df[column_order]
        pmbec_df = pmbec_df.T

        return [int(round(val * multiplier)) for val in pmbec_df.as_matrix().flatten()]

    def calculate_min_int(self):
        return min(self.as_int_list())

    def calculate_max_int(self):
        return max(self.as_int_list())
