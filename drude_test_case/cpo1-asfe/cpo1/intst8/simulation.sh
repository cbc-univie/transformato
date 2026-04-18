#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fep





path=$1

cd ${path}
pwd
hostname


for i in {1..5};
do 
istep=lig_in_vacuum
python openmm_run.py -env vacuum -odcd run_${i}/${istep}.dcd &> run_${i}/vacuum_out.log
done 



for i in {1..5};
do 
istep=lig_in_waterbox
python openmm_run.py -env waterbox -odcd run_${i}/${istep}.dcd &> run_${i}/waterbox_out.log
done 
