Transformato
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/wiederm/transformato.png)](https://travis-ci.org/wiederm/transformato)
[![codecov](https://codecov.io/gh/wiederm/transformato/branch/master/graph/badge.svg)](https://codecov.io/gh/wiederm/transformato/branch/master)

Workflow to set up a staged equilibrium sampling protocol for relative free energy calculations of ligands with a common core scaffold. The package needs to be used with output from CHARMM-GUI (www.charmm-gui.org).

### How to use transformato for binding free energy calculations

## Setting up free energy calculations
Set up your set of ligands in a waterbox and in the binding pocket of a protein using CHARMM-GUI and download the system for openMM. Put the data in the transformato/data folder. There is already an example for 2OJ9 in there, the new system should be set up the same way. E.g., if your new system is called 3BYL with ligand UN1 and UN2 you should generate two directories: transformato/data/3BYL-UN1 and transformato/data/3BYL-UN2 and move the CHARMM-GUI output for the complex in transformato/data/3BYL-UN1/complex and for the solvated ligand in transformato/data/3BYL-UN1/waterbox. Afterwards make sure that you generate a new yaml file in transformato/config. There is already an example file shown in transformato/config/2oj9-test.yaml. Most likely you only need to modify the name and the ligand three letter code of the yaml file. Then run the equilibration step of the simulation directly in transformato/data/3BYL-UN1/complex/openMM to generate the restart file for the complex and the waterbox.
## Generating equilibrium samples along the alchemical path 
After this step you should run the transformato/notebooks/example.ipynb notebook to generate the intermediate states connecting the two ligands. There are two variables you want to adapt there: the yaml file that is loaded and where you want to write the intermediate state files to. Both is further explained in the jupyter notebook.

After the notebook has generated all the intermediate states you have to run the simulations script `simulation.sh` that is located in each intermediate state directory. This script will generate equilibrium conformations and save the trajectory.
## Calculating the free energy difference using MBAR

### Maintainers

- Marcus Wieder <marcus.wieder@choderalab.org> (MSKCC)


### Copyright

Copyright (c) 2019, Marcus Wieder


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
