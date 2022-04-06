#!/bin/bash
#SBATCH -p lgpu
#SBATCH --gres=gpu


path=$1     # path in which the simulation will start
top=$2      # top file to use
script=$3   # which script is called

cd ${path}
pwd
hostname

echo 'Path: ' ${path}
echo 'top: ' ${top}
echo 'script: ' ${script}

charmm_openmm_domdec top:${top} path:${path} -i ${script}
