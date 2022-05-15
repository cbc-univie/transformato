Transformato
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/wiederm/transformato/workflows/CI/badge.svg)](https://github.com/wiederm/transformato/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/wiederm/transformato/branch/master/graph/badge.svg)](https://codecov.io/gh/wiederm/transformato/branch/master)

Transformato is a package that helps to set up relative alchemical free energy calculations of small molecules with a common core scaffold, either for solvation free energy[^1] or binding free energy estimates. The package is designed to be used with output generated by [CHARMM-GUI](https://charmm-gui.org/).

## Theory

Transformato uses a common core of two molecules to connect the physical end states of two molecules via two separate alchemical paths. The common core is defined as the maximum topology mapping between two molecular graphs, using elements to identify matching nodes.
In the example shown below the input are two physical systems (i.e. systems that are parameterized and generated without dummy atoms or unpysical additions to the topology) which are connected via their common core structure.

![SAI_Intro_v2](https://user-images.githubusercontent.com/31651017/138690737-ebc2cdcd-ee04-459e-b291-ba59b94578f8.png)

This is done by gradually turning the starting topology and parameter set in the common core molecule. For toluene and methane, the common core structure is a CH3-X, with X as a LJ particle with default vdW parameter. The alchemical transformation is performed by first turning off the electrostatics of the molecules that will become dummy atoms (for toluene this is everything that is in the green circle and for methane it is the red highlighted hydrogen). Afterwards the vdw parameters are turned off, again for the molecules in the green circle (this step is only necessary for toluene, since for methane the X hydrogen is turned into a default LJ particle). This is done via the serial atom approach (for details, see S. Boresch et. al. [^2]).


## Example systems and setup

A good starting point are the two examples of alchemical free energy calculations given in the notebooks directory.
In general, the CHARMM-GUI output for the physical end states is needed.
For the `Loeffler`-benchmark set the necessary output is shown here:  https://github.com/wiederm/transformato-systems.
To run the simulations you also need a config file, examples for these are shown in the notebooks directory.

## Installation

It is strongly recommended having anaconda (or miniconda) installed on your system which you can use for installing the dependencies (take a look in [devtools/conda-envs/test_env.yaml](https://github.com/wiederm/transformato/blob/master/devtools/conda-envs/test_env.yaml)) it creates a conda envionment called *fep*.
You can use `python setup.py install` to directly install the code of this repo into the *fep* environment.
In addition, you will need openMM or CHARMM/domdec (CHARMM/openMM) installed, depending which of the MD engines you want to use for energy/force calculations. In the *fep* environment openMM should already available.

## How to use Transformato for alchemical free energy calculations

Transformato is a python script, that generates input for MD simulations to calculated free energies between two ligands usign the CC/SAI method based on CHARMM-GUI generated folders. The CHARMM-GUI generated folders need to be structured as mentioned in point 5 later on. The **Transformato** package uses this files to create directories containing the necessary information for running an openMM or CHARMM/domdec (CHARMM/openMM) MD simulation the so called intermediate states (*intst1*,*intst2*,...). After sampling of the different intermediate states, once of the ligand solvated in water (files called *lig_in_waterbox.* and once of the protein-ligand solvated in water (corresponding files called *lig_in_complex.*), **Transformato** analyses the trajectories with the help of the python package [pymbar](https://github.com/choderalab/pymbar). A scheme of the different steps used to perform a free energy calculation with **Transformato** is shown in the folowing Figure:

<center><img src="https://user-images.githubusercontent.com/72743318/166959443-cddda959-8432-4fd9-a836-c9271b4cd2c9.svg" width="300">


For calculating relative binding free energies between two (or more ligands) one has to obtain starting structures which can be obtained from CHARMM-GUI the following way:

### Preparations in CHARMM-GUI

1. You need one pdb file containing the protein and the corresponding ligand. The ligand needs to be at the desired position (usually the binding site). In addition you need one pdf file of only the ligand and one sdf file of the ligand.
2. The pdb file needs to be uploaded to the CHARMM-GUI server ([CHARMM-GUI](https://charmm-gui.org/) use the solution builder which can be found at the input generator section after logging in to CHARMM-GUI)
3. Now go through the CHARMM-GUI processes (you may need to provide an additional sdf file of the ligand in the second step during the *PDB Manipulation options). In the *model selection* options you need to activate *Hetero* for selecting your ligand.  Select *OpenMM* or *Charmm/OpenMM* in the last step as the engine and download the final zip folder to your local device, unzip it, open the unzipped folder and rename the folder to **complex**.
4. Repeat step 3 for the same ligand without the protein (provide a pdb file containing only the ligand), name the unzipped, downloaded folder **waterbox**.
Repeat step 3 and 4 for every other Ligand.
5. Now provide all files in the following folder structure: *ligand1/waterbox*, *ligand1/complex*, *ligand2/waterbox*, *ligand2/complex* where the folders named **waterbox** and **complex** contain the downloaded CHARMM-GUI output (the needed structure can also be seen for an example in [transformato/data/cdk2-1h1q](https://github.com/wiederm/transformato/tree/master/data/cdk2-1h1q))

### Equilibration of the system

1. Before starting the simulation, you might want to equilibrate each system. Therefore, go for each system into the corresponding openmm folder (e.g. ligand1/waterbox/openmm) and run the equilibration (you can use the start of the *README* file provided by CHARMM-GUI).
2. If the equilibration was successful, you should find a file called *step4_equilibration.rst* in the complex and waterbox subfolders

### Using TRANSFORMATO for calculating a free energy between ligand 1 and ligand 2:

#### Creating intermediate state directories for the two ligands using Transformato

1. For each mutation (lig1 –> lig2) we need a `yaml` file (comparable to the one shown in [transformato/tests/config/test-28_1h1q_rbfe.yaml](https://github.com/wiederm/transformato/blob/master/transformato/tests/config/test-28_1h1q_rbfe.yaml)). You can use this file and modify it according to your needs. 
You **must ** change the following entries for each mutation:
	
	- For `name` you need to provide the folder name of the ligand where you stored the CHARMM-GUI output as **waterbox** and **complex** folder (here: *ligand1*).
	
   The following options can be kept the same for different mutations:

	-  you need the correct `tlc` description. This corresponds to the residue name of the ligand and is in most cases UNK. You can check that by looking at the *pdb* file provided by CHARMM-GUI which you can find in *ligand1/waterbox/openmm/step3_input.pdb* (4th row)
	- the number of steps for the production run (`nstep`)
	- how often trajectories are saved (`nstdcd`), and wheter you want to use constrains or not (`cons`) 
	- If you have an integration time step (`dt`) greater then 0.001 ps you might need `cons = HBonds`
	- you can select `GPU` as `True` or `False`
	- for free energy type you can select between relatvie solvation free energy `rsfe` and relative binding free energy `rbfe`
	
2. Prepare a new folder and copy the file *submit.ipynb* (an example file can be found here: [transformato/notebooks/example_toluene_methane.ipynb](https://github.com/wiederm/transformato/blob/master/notebooks/ethane-methanol-rsfe.yaml)) into it. The *submit.ipynb* is the most important python file where all commands for generating input files for the individual MD runs are generated. A few things need to be adapted:
    - The following variable needs to be modified: `configuration = load_config_yaml(config=config,input_dir=input_dir, output_dir=folder)`, for `config` you have to provide the path to the *yaml* file, for `input_dir` you need to provide the path where all your CHARMM-GUI generated folders are and for `output_dir` you have to define where the **Transformato** generated files should be stored (and where later the production run will take place!)
    - This is the standard for creating output files `perform_mutations(configuration=configuration, i=i, mutation_list=mutation_list)`, here one can additionally define `nr_of_mutation_steps_charge` (default is 5) or for ligand 1 'nr_of_mutation_steps_cc' (default is 5)


#### Running the MD simulations

The command `run_simulation(i.output_files, engine='openMM')` in the *submit.ipynb* file runs the simulation for each intermediate state one after another on your local machine. For most of the systems it might be good to use some computing cluster. If the *submit.ipynb* file is executed at the place where one wants to create the folder one can use something like the following line: 
    
```
wd = os.getcwd()
output_files = glob.glob(wd + f"/{folder}/*/*/intst*", recursive = True)

for path in sorted(output_files):
	# because path is object not string
	print(f"Start sampling for: {path}")
	exe = subprocess.Popen(["ssh", "your clustername", "sbatch", f"{str(path)}/simulation.sh", str(path)], text=True, stdout=subprocess.PIPE )
	output = exe.stdout.read()
	print(output)
```
Alternatively, as this way is quite temperamental and prone to sudden and unexplained failure, you may simply use a shell script:

```
for {i} in ./**/**/intst**/; #This assumes you start from the directory above the initial structure directory - you may change as necessary.
do cd ${i}; # to prevent issues, it is preferred to switch to the script directory
sbash simulation.sh;
cd ../../..; # switch working directory back to original
done;

```


#### Analyzing the trajectories

The trajectories can be analyzed the following way:

```
for name in (
    configuration["system"]["structure1"]["name"],
    configuration["system"]["structure2"][following],
):

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name=name,
        engine="openMM",
        max_snapshots=10000,
        num_proc=6,
        analyze_traj_with="mda",
        show_summary=True,
    )
    print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT]")
   ```
   Which gives a free energy for each ligand to the common core (``ddG_openMM``), from calculating the difference between ``ddG_openMM`` of ligand 1 and ligand 2 one obtaines the final free energy &Delta;&Delta;G(Ligand1 --> Ligand 2).
 
 
## Maintainers

- Marcus Wieder <marcus.wieder@univie.ac.at> (University of Vienna)
- Johannes Karwounopoulos <johannes.karwounopoulos@univie.ac.at> (University of Vienna)


### Copyright

Copyright (c) 2021, Marcus Wieder, Johannes Karwounopoulos, Stefan Boresch


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


[^1]:  Wieder, M., Fleck, M., Braunsfeld, B., Boresch, S., *J. Comput. Chem.* 2022, 1. [DOI](https://doi.org/10.1002/jcc.26877)
[^2]:  Boresch, S.;  Bruckner, S. *J.Comput. Chem.* 2011,32, 2449–2458 [DOI](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21829)
