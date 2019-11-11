Transformato
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/wiederm/transformato.png)](https://travis-ci.org/wiederm/transformato)
[![codecov](https://codecov.io/gh/wiederm/transformato/branch/master/graph/badge.svg)](https://codecov.io/gh/wiederm/transformato/branch/master)

Transformato is a package that helps to set up an equilibrium sampling protocol for relative free energy calculations of small molecules with a common core scaffold. The package is designed to be used with output from CHARMM-GUI (www.charmm-gui.org).

### Theory

The thermodynamic cycle used by Transformto is shown below:

![bitmap](https://user-images.githubusercontent.com/31651017/68591079-33dacb00-0490-11ea-9ca7-9c885b6f6b3f.png)

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

The common core of mol1 and mol2 is different in bonded and non-bonded parameters (mol1 is aromatic, mol2 is not aromatic). The first step is now to define a mutation route to obtain a molecule where only the common core is 


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
