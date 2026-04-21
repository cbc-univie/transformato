from transformato import load_config_yaml
from transformato.utils import postprocessing
import sys
import os
import numpy as np
import json

# comes frome the .sh file
folder = sys.argv[1]
path = sys.argv[2]
conf = sys.argv[3]

configuration = load_config_yaml(config=conf, input_dir=path, output_dir=folder)

os.chdir(folder)

print(f"#############################")
print(f"##### Structure 1 ###########")
print(f"#############################")

name1 = configuration["system"]["structure1"]["name"]

with open("configuration.json", 'w') as file:
    json.dump(configuration, file, indent=4)

struc1 = []
runs = 5 
for run in range(1, runs + 1):

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name=name1,
        engine="openMM",
        max_snapshots=10000,
        num_proc=4,
        show_summary=True,
        multiple_runs=run,
        analyze_traj_with="mda",
)
    print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT]")
    struc1.append(ddG_openMM)

print(
    f"Final free energy for {name1} is {round(np.average(struc1),2)} +- {round(np.std(struc1),2)} of the {runs} individual runs {struc1}"
)

print(
    f"f is {f_openMM}"
)
