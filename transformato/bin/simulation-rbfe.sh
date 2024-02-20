

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fep

path=$1

cd ${path}
pwd
hostname


istep=lig_in_complex
python openmm_run.py -env complex -odcd ${istep}.dcd &> complex_out.log 


istep=lig_in_waterbox
python openmm_run.py -env waterbox -odcd ${istep}.dcd  &> waterbox_out.log 
