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
import pepdata
from pepdata.amino_acid import amino_acid_letters

def get_iedb_epitopes(epitope_lengths, positive_ratio=0.6):
    df_tcell = pepdata.iedb.tcell.load_dataframe()

    # Restrict to human
    df_tcell = df_tcell[df_tcell['Host Organism Name'].fillna('').str.contains('Homo sap')]

    # Remove self
    df_tcell = df_tcell[~df_tcell['Epitope Source Organism Name'].fillna('').str.contains(
        'homo sap', case=False)]

    # Remove allergens
    for column in ['Epitope Source Molecule Name', 'In Vivo 1 Process Type',
                   'In Vivo 2 Process Type']:
        df_tcell = df_tcell[~df_tcell[column].fillna('').str.contains('allerg', case=False)]

    # Only certain lengths
    df_tcell.rename(columns={'Epitope Linear Sequence': 'iedb_epitope'}, inplace=True)
    df_tcell['epitope_length'] = df_tcell['iedb_epitope'].apply(len)
    df_tcell = df_tcell[df_tcell.epitope_length.isin(epitope_lengths)]

    # Exclude amino acid letters like B and Z that are not specific to one amino acid
    def only_amino_acid_letters(epitope):
        return all(letter in amino_acid_letters for letter in epitope)
    df_tcell = df_tcell[df_tcell.iedb_epitope.apply(only_amino_acid_letters)]

    # Calculate the T cell positive ratio, and filter by it
    df_tcell['is_tcell_positive'] = df_tcell['Qualitative Measure'].str.startswith('Positive')
    df_tcell = df_tcell[['iedb_epitope', 'epitope_length', 'is_tcell_positive']]
    def tcell_positive_ratio(bool_list):
        return sum(bool_list) / float(len(bool_list))
    df_tcell_ratio = df_tcell.groupby(['iedb_epitope', 'epitope_length']).agg(
        {'is_tcell_positive': tcell_positive_ratio,})
    df_tcell_ratio.rename(columns={'is_tcell_positive': 'positive_ratio'}, inplace=True)
    df_tcell_ratio.reset_index(inplace=True)
    df_tcell_ratio = df_tcell_ratio[df_tcell_ratio.positive_ratio >= positive_ratio][['iedb_epitope', 'epitope_length']]
    
    assert len(df_tcell_ratio.drop_duplicates()) == len(df_tcell_ratio), \
        "No duplicates should be present"
    df_tcell_ratio.reset_index(drop=True, inplace=True)
    return df_tcell_ratio
