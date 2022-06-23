#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1


path=$1
charmm=$2

cd ${path}


input=charmm_run_complex
OMP_NUM_THREADS=8 ${charmm} -i ${input}.inp > log_complex.out

input=charmm_run_waterbox
OMP_NUM_THREADS=8 ${charmm} -i ${input}.inp > log_solv.out

