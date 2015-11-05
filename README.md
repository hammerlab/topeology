[![Build Status](https://travis-ci.org/hammerlab/topeology.svg?branch=master)](https://travis-ci.org/hammerlab/topeology) [![Coverage Status](https://coveralls.io/repos/hammerlab/topeology/badge.svg?branch=master&service=github)](https://coveralls.io/github/hammerlab/topeology?branch=master)

# Topeology

Topeology compares neoepitope sequences with epitopes from [IEDB](http://www.iedb.org/).

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

## Methodology

Topeology uses Smith-Waterman alignment to align each neoepitope with each IEDB epitope of the
same length, and returns the resultant epitope-epitope scores.
