#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1


path=$1

cd ${path}
pwd
hostname

input=charmm_run_vacuum
charmm -i ${input}.inp > log_vac.out

input=charmm_run_waterbox
charmm -i ${input}.inp > log_solv.out
