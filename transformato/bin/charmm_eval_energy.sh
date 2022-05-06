#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu


path=$1     # path in which the simulation will start
top=$2      # top file to use
script=$3   # which script is called
charmm=$4   # name of the charmm executable

cd ${path}
pwd
hostname

echo 'Path: ' ${path}
echo 'top: ' ${top}
echo 'script: ' ${script}

${charmm} top:${top} path:${path} -i ${script}
