#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fep
