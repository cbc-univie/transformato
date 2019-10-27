#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -l gpu=1

. /data/shared/software/python_env/anaconda3/etc/profile.d/conda.sh
conda activate transformato

config_file=$1
working_dir=$2
structure=$4
eval_state=$4
current_state=$5

hostname
echo ${config_file}
echo ${working_dir}
echo ${structure}
echo ${eval_state}
echo ${current_state}

python  - ${config_file} ${working_dir} ${structure} ${eval_state} ${current_state} <<END

import transformato
from simtk.openmm import XmlSerializer
from simtk.openmm.app import *
import mdtraj
import json
import sys
import logging

config_file = str(sys.argv[1])
working_dir = str(sys.argv[2])
structure = str(sys.argv[3])
eval_state = int(sys.argv[4])
current_state = int(sys.argv[5])
print(system)
print(structure)
print(eval_state)
print(current_state)

conf = transformato.load_config_yaml(config=config_file, input_dir='.', output_dir=working_dir)
logging.info('State considered for energy calculations: {}'.format(eval_state))

results_for_each_env = {}
for env in ['complex', 'waterbox']:
    print(env)
    results_for_each_env[env] = transformato.calculate_energies(env, eval_state, structure, current_state=current_state, conf=conf)

json_string = json.dumps(results_for_each_env)
f = open(f"{conf['system_dir']}/results/energy_{structure}_{current_state}_{eval_state}.json", 'w+')
f.write(json_string)

END
