#!/bin/bash
#SBATCH --gres=gpu 
#SBATCH -p gpu
#SBATCH --output=analysis.log

path=$1
mol=$2

cd ${path}

current_dir=$(pwd)

echo $current_dir
echo $molecule

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate fep 
 
time python analysis.py . ./data .data/config/${mol}.yaml > ${mol}-asfe/analysis.out 
