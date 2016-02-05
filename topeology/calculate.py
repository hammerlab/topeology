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

import pandas as pd
import itertools
import imp

from .iedb_data import get_iedb_epitopes
from .scorers import CSSWLScorer, SeqAlignScorer

def get_neoepitopes(epitope_file_path, epitope_lengths):
    """
    Expected header format: sample, epitope
    """
    df_neoepitopes = pd.read_csv(epitope_file_path, dtype=object, header=0)

    # Only certain lengths
    df_neoepitopes["epitope_length"] = df_neoepitopes["epitope"].apply(len)
    df_neoepitopes = df_neoepitopes[df_neoepitopes["epitope_length"].isin(epitope_lengths)]

    df_neoepitopes.reset_index(drop=True, inplace=True)
    return df_neoepitopes

def get_joined_epitopes(epitope_file_path, epitope_lengths):
    df_neoepitopes = get_neoepitopes(epitope_file_path=epitope_file_path,
                                     epitope_lengths=epitope_lengths)
    df_iedb_epitopes = get_iedb_epitopes(epitope_lengths=epitope_lengths)
    return df_neoepitopes.merge(df_iedb_epitopes, on="epitope_length")

def calculate_similarity_from_df(df, ignore_seqalign=False):
    """
    Given a DataFrame with epitope and iedb_epitope columns, calculate
    a score for every row.
    """
    seqalign_found = False
    try:
        if not ignore_seqalign:
            imp.find_module("pmbecalign")
            seqalign_found = True
    except ImportError:
        pass

    scorer = SeqAlignScorer() if seqalign_found else CSSWLScorer()
    df["score"] = scorer.score_multiple(df, "epitope", "iedb_epitope")

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
        "sample", "epitope", "iedb_epitope", "score"]]
