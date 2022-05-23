
from calendar import c
from tabnanny import check
import numpy as np

import MDAnalysis
import yaml

# Load necessary openmm facilities
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

class Restraint():
    def __init__(
        self,
        selligand,
        selprotein,
        pdbpath,
        k,
        shape="harmonic"
        ):
        """Class representing a restraint to apply to the system

        Keywords:
        selligand,selprotein: MDAnalysis selection strings
        k: the force (=spring) constant
        pdbpath: the path to the pdbfile underlying the topology analysis
        structure: structure1 or structure2 from the yaml
        shape: 'harmonic' or 'flat-bottom' (not yet implemented): Which potential to use for the energy function"""

        self.shape=shape
        

        if self.shape not in ["harmonic","flatbottom"]:
            raise AttributeError(f"Invalid potential shape specified for restraint: {self.shape}")

        self.topology=MDAnalysis.Universe(pdbpath)
        self.g1=self.topology.select_atoms(selligand)
        self.g2=self.topology.select_atoms(selprotein)

        self.force_constant=k
        self.shape=shape
        
    def createForce(self,common_core_names):
        """Actually creates the force, after dismissing all idxs not in the common core from the molecule
        
        Keywords:
        common_core_idxs: Array - of the common core idxs"""

        
        # Only done for g1, as it is the ligand group - theres no common core for the protein
        print(f"Before CC check: {self.g1.names}")
        
        self.g1_in_cc=MDAnalysis.AtomGroup([],self.topology)
        for atom in self.g1:
            if atom.name in common_core_names:
                self.g1_in_cc+=atom
        print(f"After CC check: {self.g1_in_cc.names}")
        self.common_core_idxs=[atom.id for atom in self.g1_in_cc]
        
        # convert MDAnalysis syntax into openMM usable idxs
        self.g1_openmm=[int(id)-1 for id in self.g1_in_cc.ids]
        self.g2_openmm=[int(id)-1 for id in self.g2.ids]

        

        self.g1pos=self.g1_in_cc.center_of_mass()
        self.g2pos=self.g2.center_of_mass()

        self.initial_distance=np.linalg.norm(self.g1pos-self.g2pos)
    

        
        
        
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

        print(f"""Restraint force (centroid/bonded, shape is {self.shape}, initial distance: {self.initial_distance}):
        Group 1 (COM: {self.g1pos}): {self.g1_in_cc}
        Group 2 (COM: {self.g2pos}): {self.g2}""")

        del self.topology # delete the enormous no longer needed universe asap
        return self.force
    def get_force(self):
        return self.force
    
    def applyForce(self,system):
        """Applies the force to the openMM system"""
        system.addForce(self.force)

def get3DDistance(pos1,pos2):
    """Takes two 3d positions as array and calculates the distance between them
    
    Parameters:
    
    pos1,pos2: 3d-Array: Positions"""
    vec=pos1-pos2
    dis=np.linalg.norm(vec)
    return dis

def GenerateExtremities(configuration,pdbpath,n_extremities,sphinner=0,sphouter=5):
    """Takes the common core and generates n extremities at the furthest point
    
    Returns a selection string of the extremeties with a sphlayer selecting type C from sphinner to sphouter
    
    Keywords:
    
    configuration: the read-in restraints.yaml

    pdbpath: path to local pdb

    n_extremities: how many extremities to generate

    sphinner: Distance to start of the sphlayer, default 0

    sphouter: Distance to end of the sphlayer, default 5"""
    ligand_topology=MDAnalysis.Universe(pdbpath)
    tlc=configuration["system"]["structure"]["tlc"]
    ccs=configuration["system"]["structure"]["ccs"]
    cc_names_selection=""
    # limit ligand group to common core by way of selection string
    
    for ccname in ccs:
        cc_names_selection+=f"name {ccname} or "
    cc_names_selection=cc_names_selection[0:-4]
    print (cc_names_selection)
    ligand_group=ligand_topology.select_atoms(f"resname {tlc} and type C and ({cc_names_selection})")
    
    
    ligand_com=ligand_group.center_of_mass()
    carbon_distances={}

    if n_extremities not in [2,"2"]:
        raise NotImplementedError("Can currently only work with 2 extremities")
    for carbon in ligand_group:
        distance_from_com=get3DDistance(ligand_com,carbon.position)
        carbon_distances[carbon.name]=distance_from_com
    # Remove furthes C from group
    
    
    cs=sorted(carbon_distances,key=carbon_distances.get,reverse=True)
    furthest_carbon=ligand_group.select_atoms(f"name {cs[0]}")[0]
    print (f"Extremity 1: {furthest_carbon.name}")
    for i in cs:
        print  (f"Furthest: {i} : {carbon_distances[i]}")
    
    lg2=ligand_group.select_atoms(f"not name {furthest_carbon.name}")
    
    for carbon in lg2:
        distance_from_com=get3DDistance(furthest_carbon.position,carbon.position)
        carbon_distances[carbon.name]=distance_from_com

    cs=sorted(carbon_distances,key=carbon_distances.get,reverse=True)
    second_furthest_carbon=ligand_group.select_atoms(f"name {cs[0]}")[0]
    print (f"Extremity 2: {second_furthest_carbon.name}")
    for i in cs:
        print  (f"Furthest: {i} : {carbon_distances[i]}")

    sels=[f"name {furthest_carbon.name} or ((sphlayer {sphinner} {sphouter} name {furthest_carbon.name} and resname {tlc}) and type C)",f"(name {second_furthest_carbon.name}) or ((sphlayer {sphinner} {sphouter} name {second_furthest_carbon.name}) and resname {tlc} and type C)"]
    print (sels)
    return sels
def CreateRestraintsFromConfig(configuration,pdbpath):
    """Takes the .yaml config and returns the specified restraints
    
    Keywords:
    config: the loaded yaml config (as returned by yaml.SafeLoader() or similar"""

    tlc=configuration["system"]["structure"]["tlc"]
    
    restraints=[]
    # parse config arguments:
    restraint_command_string=configuration["simulation"]["restraints"].split()
    kval=3 #default k - value
    mode="simple"
    for arg in restraint_command_string:
        if "k=" in arg:
            kval=int(arg.split("=")[1])
        elif "extremities=" in arg:
            mode="extremities"
            n_extremities=int(arg.split("=")[1])

    if "auto" in restraint_command_string and mode=="simple":
        
        restraints.append(Restraint(f"resname {tlc} and type C" , f"(sphlayer 5 15 resname {tlc}) and name CA" , pdbpath,k=kval))
    
    elif "auto" in restraint_command_string and mode=="extremities":
        sels=GenerateExtremities(configuration,pdbpath,n_extremities)
        for selection in sels:
            restraints.append(Restraint(selection , f"(sphlayer 3 10 resname {tlc}) and name CA" , pdbpath,k=kval))
    if "manual" in restraint_command_string:
        manual_restraint_list=configuration["simulation"]["manualrestraints"].keys()
        for key in manual_restraint_list:
            restraint=configuration["simulation"]["manualrestraints"][key]
            restraints.append(Restraint(restraint["group1"],restraint["group2"],pdbpath,k=restraint["k"]))
    else:
        raise AttributeError(f"Error: demanded restraint type not recognized :{configuration['simulation']['restraints']}")

    return restraints


def write_restraints_yaml(path,system,config):
    """Takes the full config as read in in utils.py, the information for the intstate and writes the restraints.yaml
    
    path: the path to write to (e.g. ./combinedstructure/structure/intst2/restraints.yaml
    system: the system object from state.py"""
    from transformato.mutate import cc_names_struc1, cc_names_struc2
    print(cc_names_struc1)
    restraints_dict={"system":{"structure":{"tlc":f"{system.tlc}"}},"simulation":{"restraints":f"{config['simulation']['restraints']}"}}
    if "manual" in config["simulation"]["restraints"]:
        restraints_dict["simulation"]["manualrestraints"]=config["simulation"]["manualrestraints"]
    if system.structure=="structure1":
        restraints_dict["system"]["structure"]["ccs"]=cc_names_struc1
    elif system.structure=="structure2":
        restraints_dict["system"]["structure"]["ccs"]=cc_names_struc2
    
    

    output=yaml.dump(restraints_dict)
    with open(path,"w") as resyaml:
        resyaml.write(output)
        resyaml.close()
    
# Testing when running as main script
if __name__=="__main__":
    print("Initialised as main program - runing unit tests on data/2OJ9")
    

    # generate a limited config for testing purposes
    configuration={"simulation":{},"system":{"structure1":{},"structure2":{}}}
    configuration["simulation"]["restraints"]="auto"
    configuration["system"]["structure1"]["tlc"]="BMI"
    configuration["system"]["structure2"]["tlc"]="UNK"
    restraintlist=CreateRestraintsFromConfig(configuration,"../data/2OJ9-original/complex/openmm/step3_input.pdb")
    for restraint in restraintlist:
        restraint.createForce(common_core_idxs=[4824,4825,4826])
