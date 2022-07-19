"""
Manages restraints and the forces created from them.

Unless you directly want to modify functionality in here, you should not need to call here directly.
Define your restraints in the config.yaml
"""
from calendar import c
from multiprocessing.sharedctypes import Value
from tabnanny import check
import numpy as np

import MDAnalysis
import yaml
import logging

# Load necessary openmm facilities
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

logger=logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

class Restraint():
    """Class representing a restraint to apply to the system

        Intermediary step - openMM force objects are generated directly from this class and applied to the system.

        Args:
            selligand,selprotein (MDAnalysis selection string): The atoms on whose center of mass the restraint will bei applied
            k: the force (=spring) constant
            pdbpath: the path to the pdbfile underlying the topology analysis
            shape: 'harmonic' or 'flatbottom'. Defines the shape of the harmonic energy potential.
            wellsize: Defines the well-size in a flat-bottom potential. Defaults to 0.1 nanometers.
            """
    def __init__(
        self,
        selligand:str,
        selprotein:str,
        pdbpath:str,
        k:float=3,
        shape:str="harmonic",
        wellsize:float=0.05
        ):
        """Class representing a restraint to apply to the system

        Intermediary step - openMM force objects are generated directly from this class and applied to the system.

        Args:
            selligand,selprotein (MDAnalysis selection string): MDAnalysis selection strings
            k: the force (=spring) constant
            pdbpath: the path to the pdbfile underlying the topology analysis
            shape: 'harmonic' or 'flatbottom'. Defines the shape of the harmonic energy potential.
            wellsize: Defines the well-size in a flat-bottom potential. Defaults to 0.1 nanometers.
            """

        self.shape=shape
        

        if self.shape not in ["harmonic","flatbottom"]:
            raise AttributeError(f"Invalid potential shape specified for restraint: {self.shape}")

        self.topology=MDAnalysis.Universe(pdbpath)
        self.g1=self.topology.select_atoms(selligand)
        self.g2=self.topology.select_atoms(selprotein)

        self.force_constant=k
        self.shape=shape
        self.wellsize=wellsize
        
    def createForce(self,common_core_names):
        """Actually creates the force, after dismissing all idxs not in the common core from the molecule.
        
        Args:
            common_core_names (array[str]):  - Array with strings of the common core names. Usually provided by the restraint.yaml file.
            
        Returns:
            self.force: An openMM force object representing the restraint bond.
        """

        
        # Only done for g1, as it is the ligand group - theres no common core for the protein
        logger.debug(f"Before CC check: {self.g1.names}")
        
        self.g1_in_cc=MDAnalysis.AtomGroup([],self.topology)
        for atom in self.g1:
            if atom.name in common_core_names:
                self.g1_in_cc+=atom
        logger.debug(f"After CC check: {self.g1_in_cc.names}")
        self.common_core_idxs=[atom.id for atom in self.g1_in_cc]
        
        # convert MDAnalysis syntax into openMM usable idxs
        self.g1_openmm=[int(id) for id in self.g1_in_cc.ix]
        self.g2_openmm=[int(id) for id in self.g2.ix]

        logger.debug(f"G1 openMM ix: {self.g1_openmm}\nG2 openMM ix: {self.g2_openmm}")

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
                
            self.force.addBond([0,1],[self.force_constant,self.initial_distance/10])

            logger.info(f"""Restraint force (centroid/bonded, shape is {self.shape}, initial distance: {self.initial_distance}, k={self.force_constant}""")
        elif self.shape=="flatbottom":
            # create force with flat-bottom potential

            self.stepfunction=self.GenerateContinuous1DStepFunction(self.wellsize)

            self.force=CustomCentroidBondForce(2,"stepfunction((distance(g1,g2)-r0))*(distance(g1,g2)-r0)^2") # = 0 or 1 and the the harmonic potential if above the limits

            self.force.addTabulatedFunction(name="stepfunction",function=self.stepfunction)
            self.force.addPerBondParameter("k")
            self.force.addPerBondParameter("r0")
            self.force.addGroup(self.g1_openmm)
            self.force.addGroup(self.g2_openmm)
                
            self.force.addBond([0,1],[self.force_constant,self.initial_distance/10])

            logger.info(f"Force created: Flat-bottom potential. Initial distance: {self.initial_distance}. k: {self.force_constant}. Well-size: {self.wellsize}")

        else:
            raise NotImplementedError(f"Cannot create potential of shape {self.shape}")
        
        logger.debug(f"""
        Group 1 (COM: {self.g1pos}): {[self.g1_in_cc.names]}
        Group 2 (COM: {self.g2pos}): {[self.g2.names]}""")

        del self.topology # delete the enormous no longer needed universe asap
        return self.force

    def GenerateContinuous1DStepFunction(self,permitted_distance:float):
        """Creates and returns a openMM continous1DStepFunction
        
        Very simply, this creates a TabulatedFunction that takes the distance of the atom group COM to its origin, and if greater, returns 1. Otherwise, it returns 0. This function does not need to be called directly.
        It is called by the restraint generator if a flat-bottom restraint is requested.

        Args:
            

            permitted_distance: How large the well size is (meaning the 3D radius). Give in nanometers.

        Returns:
            An openMM Continuous1DFunction object representing the StepForce
        """

        # OpenMM docs: Function is assumed to be zero for x < min. Values inside range are interpolated using spline.

        
        value_table=[1,1]

        return openmm.Continuous1DFunction(value_table,min=permitted_distance,max=1500) # If there is a system larger than 1500 nm, I'm scared

        
    
    def get_force(self):
        return self.force
    
    def applyForce(self,system):
        """Applies the force to the openMM system
        
        Args:
            system (openMM.system): The openMM system to which to apply the Force
        """
        system.addForce(self.force)

def get3DDistance(pos1,pos2):
    """Takes two 3d positions as array and calculates the distance between them
    
    Args:
        pos1,pos2 (3D-Array): The positions of which to calculate the distances
        
    Returns:
        distance (float): The distance between the two positions"""
    vec=pos1-pos2
    distance=np.linalg.norm(vec)
    return distance

def GenerateExtremities(configuration,pdbpath,n_extremities,sphinner=0,sphouter=5):
    """Takes the common core and generates n extremities at the furthest point
    
        Returns a selection string of the extremities with a sphlayer selecting type C from sphinner to sphouter.
        The algorithm works as follows:

        (All of these operations only operate on carbons in the common core)
        1. Take the furthest C from the center of Mass
        2. Take the furthest C from that C
        3. Until len(extremity_cores)==n_extremities: Pick the C where the sum distance from all previous cores is greatest.
        4. Generate selection strings for all extremity cores and return
    
    Args:
        configuration (dict): the read-in restraints.yaml
        pdbpath (str): path to local pdb used as base for the restraints
        n_extremities (int): how many extremities to generate. Cannot exceed number of carbons in the ligand
        sphinner (float): Distance to start of the sphlayer, default 0
        sphouter (float): Distance to end of the sphlayer, default 5
        
    Returns:
        selection_strings (array[str]): An array of MDAnalysis selection strings, representing the selected extremities and its vicinity as defined by sphlayer
        """

    
    ligand_topology=MDAnalysis.Universe(pdbpath)
    tlc=configuration["system"]["structure"]["tlc"]
    ccs=configuration["system"]["structure"]["ccs"]
    cc_names_selection=""

    # limit ligand group to common core by way of selection string
    
    for ccname in ccs:
        cc_names_selection+=f"name {ccname} or "
    cc_names_selection=cc_names_selection[0:-4]
    logger.debug (f"Common core selection string: {cc_names_selection}")
    ligand_group=ligand_topology.select_atoms(f"resname {tlc} and type C and ({cc_names_selection})")
    
    
    ligand_com=ligand_group.center_of_mass()
    carbon_distances={}

    if int(n_extremities)<2:
        raise ValueError(f"Impossible value for generation of extremities: {n_extremities}")
    if int(n_extremities)>len(ligand_group):
        raise ValueError(f"Impossible to generate extremity restraints: too many restraints demanded ({n_extremities}) vs. carbons in common core {len(ligand_group)}")

    # Algorithm to find extremities:
    # Takes center of mass of common core. Finds furthest C from that
    # For n=2: finds furthest C from the previous C
    # For n=3: finds C, where the sum of distances to the previous furthest C is highest
    # For n=n: finds C, where the sum of distances of all previously found C is highest

    extremity_cores=[] # array containing the centers of extremity
    # Find the furthest C in the common core as starting point
    for carbon in ligand_group:
        distance_from_com=get3DDistance(ligand_com,carbon.position)
        carbon_distances[carbon.name]=distance_from_com

    # Remove furthest C from group and add it to the extremity core
    
    
    cs=sorted(carbon_distances,key=carbon_distances.get,reverse=True)
    furthest_carbon=ligand_group.select_atoms(f"name {cs[0]}")[0]
    ligand_group=ligand_group.select_atoms(f"not name {furthest_carbon.name}")
    extremity_cores.append(furthest_carbon)

    # For all other extremities - iterate over core distances, repeat process
    
    for i in range(2,int(n_extremities)+1):
        logger.debug(f"Doing Extremity number {i}")

        total_distances={}

        for carbon in ligand_group:
            total_distances[carbon.name]=0

            for core in extremity_cores:
                # Calculates the sum distance to all established cores
                total_distances[carbon.name]+=get3DDistance(core.position,carbon.position)

        # sort carbons by distances, get one with furthest distance, add to extremity_cores, remove it from ligand_group

        cs=sorted(total_distances,key=total_distances.get,reverse=True)
        resulting_carbon=ligand_group.select_atoms(f"name {cs[0]}")[0]
        ligand_group=ligand_group.select_atoms(f"not name {resulting_carbon.name}")
        extremity_cores.append(resulting_carbon)
    
    logger.info(f"Cores found for extremity_cores: {[carbon.name for carbon in extremity_cores]}")




    # Create selection strings to return
    selection_strings=[]
    for core in extremity_cores:
        selection_strings.append(f"name {core.name} or ((sphlayer {sphinner} {sphouter} name {core.name} and resname {tlc}) and type C)")

    logger.debug(f"Created extremities with selectiobns: {selection_strings}")
    return selection_strings
def CreateRestraintsFromConfig(configuration,pdbpath):
    """Takes the .yaml config and returns the specified restraints
    
    Args:
        config (dict): the loaded yaml config (as returned by yaml.safe_load() or similar
        
    Returns:
        restraints (array): An array of Restraint instances
    """

    tlc=configuration["system"]["structure"]["tlc"]
    
    restraints=[]
    # parse config arguments:
    restraint_command_string=configuration["simulation"]["restraints"].split()
    # Define default values
    restraint_args={"kval":3,"mode":"simle","shape":"harmonic","wellsize":0.1}
    for arg in restraint_command_string:
        if "k=" in arg:
            kval=int(arg.split("=")[1])
            restraint_args["k"]=kval*configuration["intst"]["scaling"]
        elif "extremities=" in arg:
            restraint_args["mode"]="extremities"
            restraint_args["n_extremities"]=int(arg.split("=")[1])
        elif "shape=" in arg:
            
            restraint_args["shape"]=str(arg.split("=")[1])
        
        elif "wellsize=" in arg:
            
            restraint_args["wellsize"]=float(arg.split("=")[1])

    if "auto" in restraint_command_string and restraint_args["mode"]=="simple":
        
        restraints.append(Restraint(f"resname {tlc} and type C" , f"(sphlayer 5 15 resname {tlc}) and name CA and protein" , pdbpath,**restraint_args))
    
    elif "auto" in restraint_command_string and restraint_args["mode"]=="extremities":
        selection_strings=GenerateExtremities(configuration,pdbpath,restraint_args["n_extremities"])
        for selection in selection_strings:
            restraints.append(Restraint(selection , f"(sphlayer 3 10 ({selection})) and name CA" , pdbpath,**restraint_args))
            
    
    if "manual" in restraint_command_string:
        manual_restraint_list=configuration["simulation"]["manualrestraints"].keys()
        for key in manual_restraint_list:
            restraint=configuration["simulation"]["manualrestraints"][key]
            restraints.append(Restraint(restraint["group1"],restraint["group2"],pdbpath,k=restraint["k"]))
    

    return restraints


def write_restraints_yaml(path,system,config,current_step):
    """Takes the full config as read in in utils.py, the information for the intstate and writes the restraints.yaml
    
    Args:
        path: the path to write to (e.g. ./combinedstructure/structure/intst2/restraints.yaml
        system: the system object from state.py
        
        """
    from transformato.mutate import cc_names_struc1, cc_names_struc2
    logger.debug(cc_names_struc1)
    restraints_dict={"intst":{},"system":{"structure":{"tlc":f"{system.tlc}"}},"simulation":{"restraints":f"{config['simulation']['restraints']}"}}

    if "scaling"  in config["simulation"]["restraints"]:
        
        if current_step==1:
            lambda_value_scaling=0
        elif current_step==2:
            lambda_value_scaling=0.25
        elif current_step==3:
            lambda_value_scaling=0.5
        elif current_step==4:
            lambda_value_scaling=0.75
        else:
            lambda_value_scaling=1
    else:
        lambda_value_scaling=1

    restraints_dict["intst"]["scaling"]=lambda_value_scaling
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
    

