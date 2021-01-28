#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1

path=$1     # path in which the simulation will start
top=$2      # top file to use
script=$3   # which script is called

cd ${path}
pwd
hostname

echo 'Path: ' ${path}
echo 'top: ' ${top}
echo 'script: ' ${script}

charmm top:${top} path:${path} -i ${script}
