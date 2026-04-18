## Test case for transformato-assisted generation of SAI-intermediate states with Drude-polarizable systems

### Main script

To setup SAI intermediate states of the molecule of interest, simply run the submit.py script, providing the residue name with the -mol option.

The test-case molecule is cyclopentanol. One would run:

`python submit.py -mol cpo1`

### Input data

Beware of hardcoded path names - the data directory provided in this test-case folder is structured such that it is guaranteed to work with transformato and may be copied/modified as a template. 

What is *required* for this workflow:

 - PSF and CRD file of the solvated system, with Drude particles
 - Configuration file(s)
 - OMM helper scripts
 - an SDF file of the molecule (no Drude particles)

#### Topology and coordinates

Found at ./data/{mol}/waterbox/openmm/step3_input.(psf|crd)

These are generated with CHARMM. Also recommended: a PDB file of the solute for the creation of the SDF.

#### Configuration

Simulation parameters are defined both in ./data/{mol}/waterbox/openmm/step5_production.psf (general) and ./data/config/{mol}.yaml.

#### OMM helper scripts

Found in the folder ./data/{mol}/waterbox/openmm/, different python modules for making OMM work with CHARMM.

#### On the SDF file

As it stands, transformato requires an SDF (./data/{mol}/{mol}/solu.sdf) for the generation of an RDKit molecule object. Importantly, the atoms in the SDF must be in the same order as the input PSF located at ./data/{mol}/waterbox/openmm/step3_input.psf for correct assignment (molecular graph nodes <> PSF entries) in the Drude case. The python script sdf_maker.py generates the required SDF from a PDB file of the solute (./data/{mol}/{mol}/solu.pdb).