# Copyright (c) 2015. Mount Sinai School of Medicine
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
from .iedb_data import get_iedb_epitopes

def similarity_score(seq_a, seq_b):
    # TODO: Add actual similarity function
    return 1.0

def get_neoepitopes(epitope_lengths, epitope_file_path):
    """
    Expected format:

    Header: sample, epitope
    """
    df_neoepitopes = pd.read_csv(epitope_file_path, dtype=object, header=0)

    # Only certain lengths
    df_neoepitopes['epitope_length'] = df_neoepitopes['epitope'].apply(len)
    df_neoepitopes = df_neoepitopes[df_neoepitopes['epitope_length'].isin(epitope_lengths)]

    df_neoepitopes.reset_index(drop=True, inplace=True)
    return df_neoepitopes

def calculate_similarity(epitope_file_path, epitope_lengths=[8, 9, 10, 11]):
    df_neoepitopes = get_neoepitopes(epitope_lengths=epitope_lengths,
                                     epitope_file_path=epitope_file_path)
    df_iedb_epitopes = get_iedb_epitopes(epitope_lengths=epitope_lengths)

    df_joined = df_neoepitopes.merge(df_iedb_epitopes, on='epitope_length')

    df_joined['score'] = df_joined.apply(
        lambda row: similarity_score(row['epitope'], row['iedb_epitope']), axis=1)
    return df_joined
