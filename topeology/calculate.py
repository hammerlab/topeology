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

import sys
import pandas as pd
import itertools
import logging
import imp

from .iedb_data import get_iedb_epitopes
from .scorers import CSSWLScorer, SeqAlignScorer

logger = logging.getLogger(__name__)
stdout_handler = logging.StreamHandler(sys.stdout)
logger.addHandler(stdout_handler)
logger.setLevel(logging.INFO)

def get_neoepitopes(epitopes, epitope_lengths):
    """
    epitopes is a file path or a DataFrame.

    Expected header format for file or DataFrame: sample, epitope
    """
    if type(epitopes) == str:
        df_neoepitopes = pd.read_csv(epitopes, dtype=object, header=0)
    elif type(epitopes) == pd.DataFrame:
        # Don't modify DataFrame that's passed in
        df_neoepitopes = epitopes.copy()
    else:
        raise ValueError('Expected str or DataFrame for epitopes argument')

    # Only certain lengths
    df_neoepitopes['epitope_length'] = df_neoepitopes['epitope'].apply(len)
    df_neoepitopes = df_neoepitopes[df_neoepitopes['epitope_length'].isin(epitope_lengths)]

    df_neoepitopes.reset_index(drop=True, inplace=True)
    return df_neoepitopes

def get_joined_epitopes(epitopes, epitope_lengths, include_hla,
                        include_organism, iedb_path, data_filters):
    df_neoepitopes = get_neoepitopes(epitopes=epitopes,
                                     epitope_lengths=epitope_lengths)
    df_iedb_epitopes = get_iedb_epitopes(epitope_lengths=epitope_lengths,
                                         include_hla=include_hla,
                                         include_organism=include_organism,
                                         iedb_path=iedb_path,
                                         data_filters=data_filters)
    return df_neoepitopes.merge(df_iedb_epitopes, on='epitope_length')

def calculate_similarity_from_df(df, ignore_seqalign=False, include_wildtype=False,
                                 include_organism=False):
    """
    Given a DataFrame with epitope and iedb_epitope columns, calculate
    a score for every row.
    """
    seqalign_found = False
    try:
        if not ignore_seqalign:
            imp.find_module('pmbecalign')
            seqalign_found = True
    except ImportError:
        pass

    if not seqalign_found:
        logger.warn('Not using seqalign for scoring')

    scorer = SeqAlignScorer() if seqalign_found else CSSWLScorer()
    df['score'] = scorer.score_multiple(df, 'epitope', 'iedb_epitope')

    if include_wildtype:
        def check_mut_wt_lengths(row):
            assert len(row['epitope']) == len(row['epitope_wt']), (
                'Mutant and wildtype epitope lengths must be equal, but saw mutant '
                '%s and wildtype %s' % (row['epitope'], row['epitope_wt']))
        df.apply(check_mut_wt_lengths, axis=1)
        df['score_wt'] = scorer.score_multiple(df, 'epitope_wt', 'iedb_epitope')

    return df

def compare(epitopes, epitope_lengths=[8, 9, 10, 11], include_wildtype=False,
            include_hla=False, include_organism=False, iedb_path=None,
            data_filters=[]):
    """
    Given a neoepitope file path or DataFrame, compare each epitope with IEDB and
    return a DataFrame with resultant scores.

    Expected header format for file or DataFrame: sample, epitope

    Output columns: sample, epitope, iedb_epitope, score
    """
    df_joined = get_joined_epitopes(epitopes=epitopes,
                                    epitope_lengths=epitope_lengths,
                                    include_hla=include_hla,
                                    include_organism=include_organism,
                                    iedb_path=iedb_path,
                                    data_filters=data_filters)
    return calculate_similarity_from_df(df_joined,
                                        include_wildtype=include_wildtype,
                                        include_organism=include_organism)
