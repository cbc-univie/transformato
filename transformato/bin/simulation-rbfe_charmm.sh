#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu


path=$1
SWITCH=$2

cd ${path}


input=charmm_run_complex
charmm_openmm_domdec -i ${input}.inp > log_complex.out

input=charmm_run_waterbox
charmm_openmm_domdec -i ${input}.inp > log_solv.out

