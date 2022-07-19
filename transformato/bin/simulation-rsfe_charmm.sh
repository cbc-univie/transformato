#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu


path=$1
SWITCH=$2

cd ${path}
pwd
hostname


run_vacuum () {
input=charmm_run_vacuum
${CHARMM} -i ${input}.inp > log_vac.out
}

run_waterbox () {
input=charmm_run_waterbox
OMP_NUM_THREADS=8 ${CHARMM} -i ${input}.inp > log_solv.out
}


case ${SWITCH} in
1)
run_vacuum
;;
2)
run_vacuum
run_waterbox
;;
esac
