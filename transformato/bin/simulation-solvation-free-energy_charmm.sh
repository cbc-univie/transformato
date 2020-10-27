#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1


#. /data/shared/software/python_env/anaconda3/etc/profile.d/conda.sh
#conda activate transformato

path=$1

cd ${path}
pwd
hostname

input=charmm_lig_in_vacuum
charmm -i ${input}.inp 

input=charmm_lig_in_waterbox
charmm -i ${input}.inp 
