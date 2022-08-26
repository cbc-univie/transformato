#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu
#SBATCH --job-name=switching # Job name
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --output=switching_%j.log # Standard output and error log


source /home/johannes/miniconda3/etc/profile.d/conda.sh
conda activate endstate

export OPENMM_PRECISION='mixed'

pwd; hostname; date


python perform_correction.py
