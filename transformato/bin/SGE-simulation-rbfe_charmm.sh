#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu
#SBATCH --exclude='n00[20-26]'


path=$1
charmm=$2

cd ${path}


input=charmm_run_complex
OMP_NUM_THREADS=8 ${charmm} -i ${input}.inp > log_complex.out

input=charmm_run_waterbox
OMP_NUM_THREADS=8 ${charmm} -i ${input}.inp > log_solv.out

