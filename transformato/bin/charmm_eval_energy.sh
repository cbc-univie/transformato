#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1

path=$1
top=$2
script=$3
switch=$4

cd ${path}
pwd
hostname

charmm top:${top} switch:${switch} -i ${script}