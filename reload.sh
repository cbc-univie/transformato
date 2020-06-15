#!/bin/bash

git checkout dev-bb
git pull origin dev-bb

conda activate transformato
conda env update -n transformato --file devtools/conda-envs/test_env.yaml

python setup.py install