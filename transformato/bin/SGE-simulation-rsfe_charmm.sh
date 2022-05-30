#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1


path=$1
SWITCH=$2

cd ${path}
pwd
hostname


input=charmm_run_vacuum
charmm_openmm_domdec -i ${input}.inp > log_vac.out

input=charmm_run_waterbox
charmm_openmm_domdec -i ${input}.inp > log_solv.out



case ${SWITCH} in
1)
run_vacuum
;;
2)
run_vacuum
run_waterbox
;;
esac
=======

