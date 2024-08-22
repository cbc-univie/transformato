path=$1

cd ${path}
pwd
hostname


istep=lig_in_vacuum
python openmm_run.py -env vacuum -odcd ${istep}.dcd -sim False &> vacuum_init_energies_out.log

# test

istep=lig_in_waterbox
python openmm_run.py -env waterbox -odcd ${istep}.dcd -sim False &> waterbox_init_energies_out.log
