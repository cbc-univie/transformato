



path=$1

cd ${path}
pwd
hostname


istep=lig_in_vacuum
python openmm_run.py -env vacuum -odcd ${istep}.dcd &> vacuum_out.log



istep=lig_in_waterbox
python openmm_run.py -env waterbox -odcd ${istep}.dcd &> waterbox_out.log
