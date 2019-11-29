#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1


. /data/shared/software/python_env/anaconda3/etc/profile.d/conda.sh
conda activate transformato

path=$1

cd ${path}
pwd
hostname

input=lig_in_complex
init=lig_in_complex
pstep=lig_in_complex
istep=lig_in_complex

python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str -orst ${istep}.rst -odcd ${istep}.dcd


input=lig_in_waterbox
init=lig_in_waterbox
pstep=lig_in_waterbox
istep=lig_in_waterbox

python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str -orst ${istep}.rst -odcd ${istep}.dcd
