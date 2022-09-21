


path=$1

source  /opt/intel/oneapi/setvars.sh

cd ${path}
pwd
hostname


input=charmm_run_vacuum
OMM_NUM_THREADS=8 charmm_c47_omm_domdecgpu -i ${input}.inp > log_vac.out

input=charmm_run_waterbox
OMM_NUM_THREADS=8 charmm_c47_omm_domdecgpu -i ${input}.inp > log_solv.out
