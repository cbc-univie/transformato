#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fep




path=$1

source  /opt/intel/oneapi/setvars.sh

cd ${path}
pwd
hostname


input=charmm_run_vacuum
OMP_NUM_THREADS=8 charmm_c47_omm_domdecgpu -i ${input}.inp > log_vac.out

input=charmm_run_waterbox
OMP_NUM_THREADS=8 charmm_c47_omm_domdecgpu -i ${input}.inp > log_solv.out
