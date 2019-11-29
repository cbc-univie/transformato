#$ -S /bin/bash
#$ -M marcus.wieder@univie.ac.at
#$ -m e
#$ -j y
#$ -p -500
#$ -o /data/shared/projects/SGE_LOG/
#$ -pe smp 3

. /data/shared/software/python_env/anaconda3/etc/profile.d/conda.sh
conda activate transformato

config_file=$1
output_dir=$2
structure=$3
conformations=$4
current_state=$5

hostname
echo ${config_file}
echo ${output_dir}
echo ${structure}
echo ${conformations}
echo ${current_state}

python  - ${config_file} ${output_dir} ${structure} ${conformations} ${current_state} <<END

import transformato
from simtk.openmm import XmlSerializer
from simtk.openmm.app import *
import mdtraj
import json
import sys, os
import logging

config_file = str(sys.argv[1])
output_dir = str(sys.argv[2])
structure = str(sys.argv[3])
conformations = int(sys.argv[4])
potential = int(sys.argv[5])
print(config_file)
print(structure)
print(conformations)
print(potential)

configuration = transformato.load_config_yaml(config=config_file, input_dir='.', output_dir=output_dir)
logging.info('Conformations at state {} are evaluated.'.format(conformations))

results_for_each_env = {}
for env in ['complex', 'waterbox']:
    print(env)
    results_for_each_env[env] = transformato.calculate_energies_with_potential_on_conf(env=env, potential=potential, conformations=conformations, structure_name=structure, configuration=configuration)

json_string = json.dumps(results_for_each_env)
os.makedirs(f"{configuration['system_dir']}/results/", exist_ok=True)
f = open(f"{configuration['system_dir']}/results/energy_{structure}_{potential}_{conformations}.json", 'w+')
f.write(json_string)

END
