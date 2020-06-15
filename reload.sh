#!/bin/bash

conda activate openmm

cd devtools/conda-envs/
conda env update -n openmm --file test_env.yaml
cd ../../

python setup.py install