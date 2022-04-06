#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu


path=$1
SWITCH=$2

cd ${path}
pwd
hostname


input=charmm_run_vacuum
charmm_openmm_domdec -i ${input}.inp > log_vac.out

input=charmm_run_waterbox
charmm_openmm_domdec -i ${input}.inp > log_solv.out

