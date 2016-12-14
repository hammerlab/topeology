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

from __future__ import print_function, absolute_import

import numpy as np
import pandas as pd
from mhctools.alleles import normalize_allele_name, compact_allele_name
import pepdata
from pepdata.amino_acid import amino_acid_letters

from .common import get_logger

logger = get_logger(__name__)

class DataFilter(object):
    def __init__(self, on, filter_col, case_sensitive=False):
        self.on = on
        self.filter_col = filter_col
        self.case_sensitive = case_sensitive

    def run_filter(self, df):
        return df[~df[self.filter_col].fillna(
            '').str.contains(self.on, case=self.case_sensitive)]

def is_valid_class_i_allele(allele):
    allele = str(allele)
    allele = allele.lower()
    if allele.startswith('hla-d') or allele.startswith('hla-e'):
        return False

    # All the valid IEDB alleles start with HLA
    if 'hla' not in allele:
        return False

    if 'undetermined' in allele:
        return False

    # Class I or class II are nonspecific
    if 'class i' in allele:
        return False

    allele = allele.upper()
    normalized_allele = normalize_allele_name(allele)
    if normalized_allele != allele:
        return False
    return True

def get_iedb_epitopes(epitope_lengths, positive_ratio=0.6, include_hla=False,
                      include_organism=False, iedb_path=None, data_filters=[]):
    """
    Note: include_hla and include_organism may increase the size of the DataFrame.
    """
    if iedb_path is not None:
        df_tcell = pd.read_csv(iedb_path,
                               skipinitialspace=True,
                               encoding='ISO-8859-1',
                               low_memory=False)
        logger.info('Using saved IEDB data at %s' % iedb_path)
    else:
        df_tcell = pepdata.iedb.tcell.load_dataframe()

    # Restrict to human
    df_tcell = df_tcell[df_tcell['Host Organism Name'].fillna('').str.contains('Homo sap')]

    # Remove self
    data_filters.insert(
        0, DataFilter(on='homo sap', filter_col='Epitope Source Organism Name'))

    # Filter out based on other attributes
    for data_filter in data_filters:
        old_len = len(df_tcell)
        df_tcell = data_filter.run_filter(df_tcell)
        new_len = len(df_tcell)
        logger.info('Filtering on %s (col %s): %d to %d elements' % (data_filter.on, data_filter.filter_col,
                                                            old_len, new_len))

    # Only certain lengths
    df_tcell.rename(columns={'Epitope Linear Sequence': 'iedb_epitope'}, inplace=True)
    df_tcell['epitope_length'] = df_tcell['iedb_epitope'].fillna('').apply(len)
    if epitope_lengths is not None:
        # pylint: disable=no-member
        # pylint gets confused by df_tcell.epitope_length here
        df_tcell = df_tcell[df_tcell.epitope_length.isin(epitope_lengths)]

    # Exclude amino acid letters like B and Z that are not specific to one amino acid
    def only_amino_acid_letters(epitope):
        return all(letter in amino_acid_letters for letter in epitope)
    # Only look at epitopes that are valid strings
    df_tcell = df_tcell[(df_tcell.iedb_epitope.apply(type) == unicode) |
                        (df_tcell.iedb_epitope.apply(type) == str)]
    df_tcell = df_tcell[df_tcell.iedb_epitope.apply(only_amino_acid_letters)]

    # Calculate the T cell positive ratio, and filter by it
    df_tcell['is_tcell_positive'] = df_tcell['Qualitative Measure'].str.startswith('Positive')

    # Include the HLA allele for any valid class I HLA alleles
    def hla_if_valid(allele):
        if not is_valid_class_i_allele(allele):
            return np.nan
        return compact_allele_name(allele)

    if include_hla:
        df_tcell['hla'] = df_tcell["MHC Allele Name"].apply(hla_if_valid)

    tcell_groupby_cols = ['iedb_epitope', 'epitope_length']
    tcell_keep_cols = ['iedb_epitope', 'epitope_length', 'is_tcell_positive']

    if include_organism:
        # pylint: disable=no-member
        # pylint gets confused by df_tcell.rename here
        df_tcell.rename(columns={'Epitope Source Organism Name': 'organism'}, inplace=True)
        tcell_keep_cols.append('organism')

    if include_hla:
        tcell_keep_cols.append('hla')

    df_tcell = df_tcell[tcell_keep_cols]
    def tcell_positive_ratio(bool_list):
        return sum(bool_list) / float(len(bool_list))
    df_tcell_ratio = df_tcell.groupby(tcell_groupby_cols).agg(
        {'is_tcell_positive': tcell_positive_ratio})
    df_tcell_ratio.rename(columns={'is_tcell_positive': 'positive_ratio'}, inplace=True)
    df_tcell_ratio.reset_index(inplace=True)
    df_tcell_ratio = df_tcell_ratio[df_tcell_ratio.positive_ratio >= positive_ratio][tcell_groupby_cols]
    assert len(df_tcell_ratio.drop_duplicates()) == len(df_tcell_ratio), \
        'No duplicates should be present'

    # Back to df_tcell; multiple rows for each HLA allele and organism, if applicable
    df_tcell = df_tcell.merge(df_tcell_ratio, on=tcell_groupby_cols)

    # Now that we've filtered to only epitopes with a certain T-cell positive ratio, don't include any rows
    # for HLA alleles and organisms for which there wasn't a T-cell positive result
    df_tcell = df_tcell[df_tcell.is_tcell_positive]
    df_tcell = df_tcell.drop(labels='is_tcell_positive', axis=1)

    # Might have e.g. duplicate HLA values
    df_tcell.drop_duplicates(inplace=True)
    df_tcell.reset_index(drop=True, inplace=True)
    return df_tcell
