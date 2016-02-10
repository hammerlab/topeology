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
Generate test data of a particular length.

%(prog)s --length 10 --output test.csv --epitope-lengths 8 9 10 11
"""
from __future__ import print_function, absolute_import
import argparse
from topeology import compare
from numpy.random import randint, choice
import pandas as pd

AMINO_ACID_LETTERS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                      'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                      'T', 'W', 'Y', 'V']

def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        '--length',
        nargs=1,
        type=int,
        required=True,
        help='Number of data rows to generate')
    parser.add_argument(
        '--output',
        nargs=1,
        required=True,
        help='Output file path')
    parser.add_argument(
        '--epitope-lengths',
        nargs='+',
        type=int,
        default=[8, 9, 10, 11],
        help='Allowed epitope lengths to be compared with')
    args = parser.parse_args()
    length = args.length[0]
    sample_list = randint(100, size=length)
    epitope_list = [randepitope(args.epitope_lengths) for i in range(length)]
    df = pd.DataFrame({'sample': sample_list, 'epitope': epitope_list})
    df.set_index('sample', inplace=True)
    df.to_csv(args.output[0])

def randepitope(sizes):
    ints = randint(len(AMINO_ACID_LETTERS), size=choice(sizes))
    chars = [AMINO_ACID_LETTERS[i] for i in ints]
    return ''.join(chars)

if __name__ == '__main__':
    run()
