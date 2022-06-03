#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu
path=$1
SWITCH=$2

cd ${path}
pwd
hostname


run_complex () {
input=charmm_run_complex
OMP_NUM_THREADS=8 ${CHARMM} -i ${input}.inp > log_complex.out
}

run_waterbox () {
input=charmm_run_waterbox
OMP_NUM_THREADS=8 ${CHARMM} -i ${input}.inp > log_solv.out
}


case ${SWITCH} in
1)
run_complex
;;
2)
run_complex
run_waterbox
;;
esac
