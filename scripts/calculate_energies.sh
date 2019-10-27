#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1

. /data/shared/software/python_env/anaconda3/etc/profile.d/conda.sh
conda activate sai-binding

system=$1
structure=$2
eval_state=$3
current_state=$4

hostname
echo ${system}
echo ${structure}
echo ${eval_state}
echo ${current_state}

python  - ${system} ${structure} ${eval_state} ${current_state} <<END

import sai
from simtk.openmm import XmlSerializer
from simtk.openmm.app import *
import mdtraj
import json
import sys
import logging

system = str(sys.argv[1])
structure = str(sys.argv[2])
eval_state = int(sys.argv[3])
current_state = int(sys.argv[4])
print(system)
print(structure)
print(eval_state)
print(current_state)

conf = transformato.load_config_yaml(system)
logging.info('State considered for energy calculations: {}'.format(eval_state))

results_for_each_env = {}
for env in ['complex', 'waterbox']:
    print(env)
    results_for_each_env[env] = transformato.calculate_energies(env, eval_state, structure, current_state=current_state, conf=conf)

json_string = json.dumps(results_for_each_env)
f = open(f"{conf['system_dir']}/results/energy_{structure}_{current_state}_{eval_state}.json", 'w+')
f.write(json_string)

END
