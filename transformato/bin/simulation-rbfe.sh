

source ~/anaconda3/etc/profile.d/conda.sh
conda activate fep

path=$1

cd ${path}
pwd
hostname

input=lig_in_complex
init=lig_in_complex
pstep=lig_in_complex
istep=lig_in_complex
irst=lig_in_complex
orst=lig_in_complex_rst
python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${irst}.rst -orst ${irst} -odcd ${istep}.dcd &> complex_out.log 

input=lig_in_waterbox
init=lig_in_waterbox
pstep=lig_in_waterbox
istep=lig_in_waterbox
irst=lig_in_waterbox
orst=lig_in_waterbox_rst
python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${irst}.rst -orst ${irst} -odcd ${istep}.dcd  &> waterbox_out.log 
