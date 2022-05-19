from xml.dom.minidom import Attr
import numpy as np

import MDAnalysis

# Load necessary openmm facilities
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

class Restraint():
    def __init__(self,g1,g2,k=3,shape="harmonic"):
        """Class representing a restraint to apply to the system

        Keywords:
        g1,g2: An array of idxs representing the groups to be bonded
        k: the force (=spring) constant
        shape: 'harmonic' or 'flat-bottom' (not yet implemented): Which potential to use for the energy function"""

        self.g1=g1
        self.g2=g2
        self.force_constant=k
        self.shape=shape
        if self.shape not in ["harmonic","flatbottom"]:
            raise AttributeError(f"Invalid potential shape specified for restraint: {self.shape}")


def CreateRestraintsFromConfig(configuration,pdbpath):
    """Takes the .yaml config and returns the specified restraints
    
    Keywords:
    config: the loaded yaml config (as returned by yaml.SafeLoader() or similar"""

    topology=MDAnalysis.Universe(pdbpath)

    if configuration["simulation"]["restraints"]=="auto":
        CreateSimpleRestraints(configuration,topology,"structure1")
    
    else:
        raise AttributeError(f"Error: demanded restraint type not supported :{configuration['simulation']['restraints']}")


def CreateSimpleRestraints(configuration,topology,structure):

    tlc=configuration["system"][structure]["tlc"]
    g1=topology.select_atoms(f"resname {tlc} and name C**")
    g2=topology.select_atoms(f"(sphlayer 5 15 resname {tlc}) and name CA") # nearby carbons for ligand

    g1_pos=g1.center_of_mass()
    g2_pos=g2.center_of_mass()

    initial_distance=np.linalg.norm(g1_pos-g2_pos)
    print(f"""Restraint (centroid/bonded, initial distance: {initial_distance}) created:
    Group 1 (COM: {g1_pos}): {g1}
    Group 2 (COM: {g2_pos}): {g2}""")

# Testing when running as main script
if __name__=="__main__":
    print("Initialised as main program - runing unit tests")
    testrestraint=Restraint(22,24,5,"harmonic")

    # generate a limited config for testing purposes
    configuration={"simulation":{},"system":{"structure1":{},"structure2":{}}}
    configuration["simulation"]["restraints"]="auto"
    configuration["system"]["structure1"]["tlc"]="BMI"
    configuration["system"]["structure2"]["tlc"]="UNK"
    CreateRestraintsFromConfig(configuration,"../data/2OJ9-original/complex/openmm/step3_input.pdb")

