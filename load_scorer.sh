#!/bin/bash

set -o errexit

source ./ENV.sh

python ./load_scorer.py
