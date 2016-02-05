[![Build Status](https://travis-ci.org/hammerlab/topeology.svg?branch=master)](https://travis-ci.org/hammerlab/topeology) [![Coverage Status](https://coveralls.io/repos/hammerlab/topeology/badge.svg?branch=master&service=github)](https://coveralls.io/github/hammerlab/topeology?branch=master)

# Topeology

Topeology compares neoepitope sequences with epitopes from [IEDB](http://www.iedb.org/).

This software is **not** currently ready for production use.

## Example

From the command line:

```sh
topeology --input epitopes.csv --epitope-lengths 8 9 10 11 > scores.csv
```

In Python:

```python
from topeology import compare
output_dataframe = compare('epitopes.csv')
```

Input looks like:

| sample      | epitope
| ------      | -------
| 001         | AAALPGKCGV

Output looks like:

| sample      | epitope        | iedb_epitope    | score
| ------      | -------        | ------------    | -----
| 001         | AAALPGKCGV     | EFKEFAAGRR      | 2.38

## Installation

You can install topeology using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```sh
pip install topeology
```

Currently, topeology wraps [seq-align] in a C extension, which will be installed if [seq-align]
is installed. (Otherwise, topeology reverts to using another scorer.)

To set up [seq-align] for topeology: follow [seq-align]'s installation instructions, and then set
`SEQ_ALIGN_PATH` to your [seq-align] installation directory before installing topeology.

If topeology is already installed with [seq-align], perform the above steps and then run
`pip install topeology --upgrade --no-deps --force-reinstall`.

## Methodology

Topeology uses Smith-Waterman alignment to align each neoepitope with each IEDB epitope of the
same length, and returns the resultant epitope-epitope scores. Only position 3 to the penultimate
amino acid are considered.

This software uses the following libraries for Smith-Waterman alignment:

- [seq-align]
- [Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

[seq-align]: https://github.com/noporpoise/seq-align