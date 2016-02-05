# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Compare neoepitopes with epitopes from IEDB.

%(prog)s --input epitopes.csv --epitope-lengths 8 9 10 > scores.csv
"""
from __future__ import print_function, absolute_import
import argparse
from topeology import compare

def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        '--input',
        nargs=1,
        required=True,
        help='CSV file containing predicted epitopes and corresponding HLA alleles')
    parser.add_argument(
        '--epitope-lengths',
        nargs='+',
        type=int,
        default=[8, 9, 10, 11],
        help='CSV file containing predicted epitopes and corresponding HLA alleles')
    args = parser.parse_args()
    if len(args.input) > 1:
        raise ValueError('Only a single --input file is allowed.')
    compare_df = compare(epitope_file_path=args.input[0],
                         epitope_lengths=args.epitope_lengths)
    print(compare_df.to_csv(index=False))
    
if __name__ == '__main__':
    run()
