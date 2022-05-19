
import numpy as np

import MDAnalysis

# Load necessary openmm facilities
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

class Restraint():
    def __init__(self,sel1,sel2,pdbpath,k=3,shape="harmonic"):
        """Class representing a restraint to apply to the system

        Keywords:
        g1,g2: An MDAnalysis selection string
        k: the force (=spring) constant
        shape: 'harmonic' or 'flat-bottom' (not yet implemented): Which potential to use for the energy function"""

        self.shape=shape

        if self.shape not in ["harmonic","flatbottom"]:
            raise AttributeError(f"Invalid potential shape specified for restraint: {self.shape}")

        self.topology=MDAnalysis.Universe(pdbpath)
        self.g1=self.topology.select_atoms(sel1)
        self.g2=self.topology.select_atoms(sel2)

        self.g1pos=self.g1.center_of_mass()
        self.g2pos=self.g2.center_of_mass()

        self.initial_distance=np.linalg.norm(self.g1pos-self.g2pos)
    

        # convert MDAnalysis syntax into openMM usable idxs
        self.g1_openmm=[int(id)-1 for id in self.g1.ids]
        self.g2_openmm=[int(id)-1 for id in self.g2.ids]

        self.force_constant=k
        self.shape=shape
        
        if self.shape=="harmonic":
            # create force with harmonic potential
            self.force=CustomCentroidBondForce(2,"0.5*k*(distance(g1,g2)-r0)^2")
            self.force.addPerBondParameter("k")
            self.force.addPerBondParameter("r0")
            self.force.addGroup(self.g1_openmm)
            self.force.addGroup(self.g2_openmm)
                
            self.force.addBond([0,1],[self.force_constant,self.initial_distance])

            
        elif self.shape=="flatbottom":
            raise NotImplementedError("Cant create flatbottom potential, as it has not yet been implemented")

        print(f"""Restraint force (centroid/bonded, shape is {self.shape}, initial distance: {self.initial_distance}) created:
        Group 1 (COM: {self.g1pos}): {self.g1_openmm}
        Group 2 (COM: {self.g2pos}): {self.g2_openmm}""")
    
    def get_force(self):
        return self.force
    
    
def CreateRestraintsFromConfig(configuration,pdbpath):
    """Takes the .yaml config and returns the specified restraints
    
    Keywords:
    config: the loaded yaml config (as returned by yaml.SafeLoader() or similar"""

    tlc=configuration["system"]["structure1"]["tlc"]
    restraints=[]
    if configuration["simulation"]["restraints"]=="auto":
        
        restraints.append(Restraint(f"resname {tlc} and name C**" , f"(sphlayer 5 15 resname {tlc}) and name CA" , pdbpath))
    
    else:
        raise AttributeError(f"Error: demanded restraint type not supported :{configuration['simulation']['restraints']}")





# Testing when running as main script
if __name__=="__main__":
    print("Initialised as main program - runing unit tests")
    

    # generate a limited config for testing purposes
    configuration={"simulation":{},"system":{"structure1":{},"structure2":{}}}
    configuration["simulation"]["restraints"]="auto"
    configuration["system"]["structure1"]["tlc"]="BMI"
    configuration["system"]["structure2"]["tlc"]="UNK"
    CreateRestraintsFromConfig(configuration,"../data/2OJ9-original/complex/openmm/step3_input.pdb")

