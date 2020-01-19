Transformato
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/wiederm/transformato.png)](https://travis-ci.org/wiederm/transformato)
[![codecov](https://codecov.io/gh/wiederm/transformato/branch/master/graph/badge.svg)](https://codecov.io/gh/wiederm/transformato/branch/master)

Transformato is a package that helps to set up an equilibrium sampling protocol for relative free energy calculations of small molecules with a common core scaffold. The package is designed to be used with output from CHARMM-GUI (www.charmm-gui.org).

### Theory

The thermodynamic cycle used by Transformto is shown below:

![drawing-1](https://user-images.githubusercontent.com/31651017/72672192-dc8c4480-3a56-11ea-9906-67bd27810c9e.png)

The strategy is to turn atoms of two molecules to dummy atoms (shown in green) to obtain a common core (cc) between the two molecules (in this case a CH3 group) and subsequently mutate the common core of molecule 1 (cc1) to the common core of molecule 2 (cc2) to obtain equivalen states. This closes the thermodynamic cycle. 
In the depicted example the additional step from cc1 to cc2 might not be necessary since both cc1-CH3 and cc2-CH3 might be equivalent in their bonded and non-bonded terms, but this is a special case. For more general mutations the change in bonded and non-bonded terms will be necessary. The mutation strategy used by Transformato can deal with changes in bonded and non-bonded parameters of the common core.   


A more in detail description of the mutation strategy used by Transformto will be described using the the following two molecules:
![bitmap](https://user-images.githubusercontent.com/31651017/68553941-6e5b4e00-0425-11ea-8b81-065e013276c2.png)

The general approach is:
- find the common maximum (connected) substructure (CMS) based on elements (ignoring bond and atom types)
  - this defines the commen core (cc) for both small molecules resulting in cc1 and cc2. The commen core can be different in atom types, charges and bonded parameters but has a mapping between the atoms.
- mutate the atoms that are not in the CMS.
  - charges are linearly scaled to zero in multiple lambda states.
  - LJ parameters are set to zero using a mutation strategy called serial atom insertion (SAI). This leads to either LJ fully intact or completly turned off for a single atom outside of the common core. At each LJ lambda protocol step the LJ terms of at least one atom are turned off.
- linearly scaling of the bonded and non-bonded parameters of the common core.

Using these transformations we obtain the same common core between two different molecules with a common substructure.


## Common maximum substructure

For our example the maximum common substructure is marked in read on both molecules.
![bitmap](https://user-images.githubusercontent.com/31651017/68554631-c2683180-0429-11ea-8846-ce303fa933b9.png)

The common core of mol1 and mol2 is different in bonded and non-bonded parameters (mol1 is aromatic, mol2 is not aromatic). The first step is now to define a mutation route to obtain a molecule where only the common core is interacting with the environment.


### How to use transformato for binding free energy calculations

## Setting up free energy calculations

For 'endpoints' are needed to run a protein binding free energy simulation. 
1. Ligand 1 in a waterbox
2. Ligand 2 in a waterbox
3. Ligand 1 in complex
4. Ligand 2 in complex

An example with the ligands shown above is deposited in transformato/data.
The directory with 2OJ9-e1 has the enol form of the ligand and 2OJ9-e2 the keto form of the ligand.

All these directories are directly taken from CHARMM-GUI. After you solvate either the ligand or the ligand+protein and download the system for openMM you have to run a short minimization using the scripts obtaine from CHARMM-GUI.

To use the Transformato free energy set up the next step is to set up a config file for you run.
The config file for 2OJ9-e1 and 2OJ9-e2 is shown in transformato/config/2oj9_example.yaml.
What is important here is that the tree letter code (tlc) used in the psf is the same as in the config file.


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
