"""
Manages restraints and the forces created from them.

Basic function:
-------------------

From the 'restraints' section in the config.yaml, the intermediate state factory creates an appropriate restraints.yaml in the
intermediate state directories.

When the simulation is run, openMM_run.py checks for the existence of these restraints.yaml. If it finds one, a number of Restraint() instances commensurate
to the number of desired restraints is created. These each create an appropriate openMM force as defined by their parameters. That force ist then applied to the openMM simulation.

Usage:
--------

Unless you directly want to modify functionality, you should not need to call here directly.
Define your restraints in the config.yaml.

.. hint::
    For further usage instructions, see the 'System Setup' section of the main documentation.

"""

import logging

import MDAnalysis
import numpy as np
import yaml

# Load necessary openmm facilities
from openmm import CustomCentroidBondForce
from openmm.unit import angstrom

logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)


class Restraint:
    def __init__(
        self,
        selligand: str,
        selprotein: str,
        topology: MDAnalysis.Universe,
        k: float = 3,
        shape: str = "harmonic",
        wellsize: float = 0.05,
        **kwargs,
    ):
        """Class representing a restraint to apply to the system.

        Intermediary step - openMM force objects are generated directly from this class and applied to the system.

        Raises:
            AttributeError: If the potential shape requested is not implemented.

        Args:
            selligand,selprotein (MDAnalysis selection string): MDAnalysis selection strings
            k: the force (=spring) constant applied to the potential energy formula. See the 'System Setup' section for details.
            topology: the MDAnalysis universe used to generate restraint geometries
            shape: one of 'harmonic', 'flatbottom', 'flatbottom-oneside-sharp' or 'flatbottom-twoside'. Defines the shape of the harmonic energy potential.
            wellsize: Defines the well-size in a two-sided flat-bottom potential. Defaults to 0.05 nanometers.
            kwargs: Catcher for additional restraint_args

        """

        self.shape = shape

        if self.shape not in [
            "harmonic",
            "flatbottom",
            "flatbottom-oneside-sharp",
            "flatbottom-oneside",
            "flatbottom-twoside",
        ]:
            raise AttributeError(
                f"Invalid potential shape specified for restraint: {self.shape}"
            )

        self.topology = topology
        self.g1 = self.topology.select_atoms(selligand)
        self.g2 = self.topology.select_atoms(selprotein)

        self.force_constant = k
        self.shape = shape
        self.wellsize = wellsize * angstrom
        self.kwargs = kwargs

    def createForce(self, common_core_names):
        """Actually creates the force, after dismissing all idxs not in the common core from the molecule.

        Args:
            common_core_names (array[str]):  - Array with strings of the common core names. Usually provided by the restraint.yaml file.

        Returns:
            openMM.CustomCentroidBondForce: An openMM force object representing the restraint bond.
        """

        # Only done for g1, as it is the ligand group - theres no common core for the protein
        logger.debug(f"Before CC check: {self.g1.names}")

        self.g1_in_cc = MDAnalysis.AtomGroup([], self.topology)
        for atom in self.g1:
            if atom.name in common_core_names:
                self.g1_in_cc += atom
        logger.debug(f"After CC check: {self.g1_in_cc.names}")
        self.common_core_idxs = [atom.id for atom in self.g1_in_cc]

        # convert MDAnalysis syntax into openMM usable idxs
        self.g1_openmm = [int(id) for id in self.g1_in_cc.ix]
        self.g2_openmm = [int(id) for id in self.g2.ix]

        logger.debug(f"G1 openMM ix: {self.g1_openmm}\nG2 openMM ix: {self.g2_openmm}")

        self.g1pos = self.g1_in_cc.center_of_mass()
        self.g2pos = self.g2.center_of_mass()

        self.initial_distance = np.linalg.norm(self.g1pos - self.g2pos)

        if (
            not "r0" in self.kwargs.keys()
        ):  # default value - equilibrium distance is the initial distance. Otherwise, overriden with the r0 specified
            self.r0 = self.initial_distance * angstrom
        else:
            self.r0 = self.kwargs["r0"] * angstrom

        if self.shape == "harmonic":
            # create force with harmonic potential
            self.force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
            self.force.addPerBondParameter("k")
            self.force.addPerBondParameter("r0")
            self.force.addGroup(self.g1_openmm)
            self.force.addGroup(self.g2_openmm)

            self.force.addBond([0, 1], [self.force_constant, self.r0])

            logger.info(
                f"""Restraint force (centroid/bonded, shape is {self.shape}, initial distance: {self.initial_distance}, k={self.force_constant}"""
            )
        elif self.shape in ["flatbottom", "flatbottom-oneside"]:
            # create force with flat-bottom potential identical to the one used in openmmtools
            logger.debug("creating flatbottom-1side-soft")
            self.force = CustomCentroidBondForce(
                2, "step(distance(g1,g2)-r0) * (k/2)*(distance(g1,g2)-r0)^2"
            )

            self._add_flatbottom_parameters()

        elif self.shape == "flatbottom-oneside-sharp":
            # create force with flat-bottom potential
            logger.debug("creating flatbottom-1side-sharp")
            self.force = CustomCentroidBondForce(
                2, "step(distance(g1,g2)-r0) * (k/2)*(distance(g1,g2))^2"
            )

            self._add_flatbottom_parameters()

        elif self.shape == "flatbottom-twoside":
            # create force with flat-bottom potential
            logger.debug("creating flatbottom-2side-sharp")
            self.force = CustomCentroidBondForce(
                2, "step(abs(distance(g1,g2)-r0)-w)*k*(distance(g1,g2)-r0)^2"
            )

            self._add_flatbottom_parameters()

        else:
            raise NotImplementedError(f"Cannot create potential of shape {self.shape}")

        logger.debug(
            f"""
        Group 1 (COM: {self.g1pos}): {[self.g1_in_cc.names]}
        Group 2 (COM: {self.g2pos}): {[self.g2.names]}"""
        )

        del self.topology  # delete the enormous no longer needed universe asap
        return self.force

    def get_force(self):
        return self.force

    def _add_flatbottom_parameters(self):
        """Takes care of the initial setup of an openMM CustomCentroidBondForce with a flatbottom shape

        Args:
            force: An openMM CustomCentroidBondForce object
        """

        self.force.addPerBondParameter("k")
        self.force.addPerBondParameter("r0")

        self.force.addGroup(self.g1_openmm)
        self.force.addGroup(self.g2_openmm)

        if (
            "twoside" in self.shape
        ):  # two - sided flatbottoms gain an additional parameter for the wellsize
            self.force.addPerBondParameter("w")
            self.force.addBond([0, 1], [self.force_constant, self.r0, self.wellsize])
        else:
            self.force.addBond([0, 1], [self.force_constant, self.r0])

        logger.info(
            f"Restraint force (centroid/bonded), shape is {self.shape}, Parameters: {self.force.getBondParameters(0)} "
        )

    def applyForce(self, system):
        """Applies the force to the openMM system

        Args:
            system (openMM.system): The openMM system to which to apply the Force
        """
        system.addForce(self.force)


def get3DDistance(pos1, pos2):
    """Takes two 3d positions as array and calculates the distance between them

    Args:
        pos1,pos2 (3D-Array): The positions of which to calculate the distances

    Returns:
        float: The distance between the two positions
    """
    vec = pos1 - pos2
    distance = np.linalg.norm(vec)
    return distance


def generate_simple_selection(configuration):
    """Takes the common core and selects surrounding carbon-alphas

    This ensures that the initial simple restraints on both sides is identical

    Args:
        configuration (dict): the read-in restraints.yaml
        

    Returns:
        str: An MDAnalysis selection string, representing the carbon-alphas surrounding the cores.
    """

    tlc = configuration["system"]["structure"]["tlc"]
    ccs = configuration["system"]["structure"]["ccs"]
    cc_names_selection = ""

    # limit ligand group to common core by way of selection string

    for ccname in ccs:
        cc_names_selection += f"name {ccname} or "
    cc_names_selection = cc_names_selection[0:-4]

    selstr = f"(sphlayer 5 15 ({cc_names_selection})) and name CA and protein"
    logger.debug(f"Common core selection string: {selstr}")
    return selstr


def generate_extremities(configuration, topology, n_extremities, sphinner=0, sphouter=5):
    """Takes the common core and generates n extremities at the furthest point

        Returns a selection string of the extremities with a sphlayer (see MDAnalysis docs) selecting type C from sphinner to sphouter.
        The algorithm works as follows:

        (**All of these operations only operate on carbons in the common core**)


        1. Take the furthest C from the center of Mass.

        2. Take the furthest C from that C.

        3. Until len(extremity_cores)==n_extremities:

            Pick the C where the sum distance from all previous cores is greatest.

            Add it to the extremity_cores.

        4. Generate selection strings for all extremity cores and return.

    Args:
        configuration (dict): the read-in restraints.yaml
        topology (MDAnalysis.Universe): MDAnalysis.Universe object used to generate ligand geometries
        n_extremities (int): how many extremities to generate. Cannot exceed number of carbons in the ligand
        sphinner (float): Distance to start of the sphlayer, default 0
        sphouter (float): Distance to end of the sphlayer, default 5

    Raises:
        ValueError: If an invalid amount of extremities is specified.

    Returns:
        [selection_strings,extremity_cores]: A nested array of MDAnalysis selection strings representing the selected extremities and its vicinity as defined by sphlayer and an array of the found extremity cores
    """

    ligand_topology = topology
    tlc = configuration["system"]["structure"]["tlc"]
    ccs = configuration["system"]["structure"]["ccs"]
    cc_names_selection = ""

    # limit ligand group to common core by way of selection string

    for ccname in ccs:
        cc_names_selection += f"name {ccname} or "
    cc_names_selection = cc_names_selection[0:-4]
    logger.debug(f"Common core selection string: {cc_names_selection}")
    ligand_group = ligand_topology.select_atoms(
        f"resname {tlc} and type C and ({cc_names_selection})"
    )

    ligand_com = ligand_group.center_of_mass()
    carbon_distances = {}

    if int(n_extremities) < 2:
        raise ValueError(
            f"Impossible value for generation of extremities: {n_extremities}"
        )
    if int(n_extremities) > len(ligand_group):
        raise ValueError(
            f"Impossible to generate extremity restraints: too many restraints demanded ({n_extremities}) vs. carbons in common core {len(ligand_group)}"
        )

    # Algorithm to find extremities:
    # Takes center of mass of common core. Finds furthest C from that
    # For n=2: finds furthest C from the previous C
    # For n=3: finds C, where the sum of distances to the previous furthest C is highest
    # For n=n: finds C, where the sum of distances of all previously found C is highest

    extremity_cores = []  # array containing the centers of extremity
    # Find the furthest C in the common core as starting point
    for carbon in ligand_group:
        distance_from_com = get3DDistance(ligand_com, carbon.position)
        carbon_distances[carbon.name] = distance_from_com

    # Remove furthest C from group and add it to the extremity core

    cs = sorted(carbon_distances, key=carbon_distances.get, reverse=True)
    furthest_carbon = ligand_group.select_atoms(f"name {cs[0]}")[0]
    ligand_group = ligand_group.select_atoms(f"not name {furthest_carbon.name}")
    extremity_cores.append(furthest_carbon)

    # For all other extremities - iterate over core distances, repeat process

    for i in range(2, int(n_extremities) + 1):
        logger.debug(f"Doing Extremity number {i}")

        total_distances = {}

        for carbon in ligand_group:
            total_distances[carbon.name] = 0

            for core in extremity_cores:
                # Calculates the sum distance to all established cores
                total_distances[carbon.name] += get3DDistance(
                    core.position, carbon.position
                )

        # sort carbons by distances, get one with furthest distance, add to extremity_cores, remove it from ligand_group

        cs = sorted(total_distances, key=total_distances.get, reverse=True)
        resulting_carbon = ligand_group.select_atoms(f"name {cs[0]}")[0]
        ligand_group = ligand_group.select_atoms(f"not name {resulting_carbon.name}")
        extremity_cores.append(resulting_carbon)

    logger.info(
        f"Cores found for extremity_cores: {[carbon.name for carbon in extremity_cores]}"
    )

    # Create selection strings to return
    selection_strings = []
    for core in extremity_cores:
        selection_strings.append(
            f"name {core.name} or ((sphlayer {sphinner} {sphouter} name {core.name} and resname {tlc}) and type C)"
        )

    logger.debug(f"Created extremities with selections: {selection_strings}")
    return selection_strings,extremity_cores


def create_restraints_from_config(configuration, pdbpath):
    """Takes the restraints.yaml config and returns the specified restraints.

    This is the function that should be invoked by openmm_run.py.

    Args:
        config (dict): the loaded yaml config (as returned by yaml.safe_load() or similar

    Returns:
        array: An array of Restraint instances
    """
    universe=MDAnalysis.Universe(pdbpath)

    tlc = configuration["system"]["structure"]["tlc"]

    restraints = []
    # parse config arguments and pass to restraint generation:
    restraint_command_string = configuration["simulation"]["restraints"].split()
    restraint_args = {"mode": "simple"}

    for arg in restraint_command_string:
        if "k=" in arg:
            kval = int(arg.split("=")[1])
            restraint_args["k"] = kval * configuration["intst"]["scaling"]
        elif "extremities=" in arg:
            restraint_args["mode"] = "extremities"
            restraint_args["n_extremities"] = int(arg.split("=")[1])
        elif "shape=" in arg:

            restraint_args["shape"] = str(arg.split("=")[1])

        elif "wellsize=" in arg:

            restraint_args["wellsize"] = float(arg.split("=")[1])

    if "auto" in restraint_command_string and restraint_args["mode"] == "simple":
        logger.debug("generating simple selection")
        selstr = generate_simple_selection(configuration)
        restraints.append(
            Restraint(f"resname {tlc} and type C", selstr, universe, **restraint_args)
        )

    elif "auto" in restraint_command_string and restraint_args["mode"] == "extremities":
        logger.debug("generating extremity selections")
        selection_strings, extremity_cores = generate_extremities(
            configuration, universe, restraint_args["n_extremities"]
        )
        for i,selection in enumerate(selection_strings):
            restraints.append(
                Restraint(
                    selection,
                    f"(sphlayer 3 10 ({selection})) and name CA",
                    universe,
                    ex_core=extremity_cores[i],
                    **restraint_args
                )
            )
        
        # At this point only automatic ex-restraints exist

        ex_cores=[restraint.kwargs["ex_core"] for restraint in restraints]

        logger.info(15*"-"+f"Available cores for assignment: {len(ex_cores)}")
        

        def assign_atom(atom:MDAnalysis.core.groups.Atom):
            """
            Inner function to assign an atom involved in multiple restraint ligand groups to a single one.
            
            The restraint with the geometrically closest core will be chosen and the atom deleted from all others.

            
            Args:
            atom: The duplicate atom to reassign.
            """

            distances=dict()
            

            for core in ex_cores:
                distances[core]=get3DDistance(atom.position,core.position)
                
                sorted_distances=dict(sorted(distances.items(),key=lambda x:x[1]))

            logger.info(f"Sorted Distances: {sorted_distances}")
            closest=list(sorted_distances.keys())[0]
            for restraint in restraints:
                if restraint.kwargs["ex_core"].ix!=closest.ix:
                    logger.info(f"Removing {atom}")
                    logger.info(f"Old: {restraint.g1.names}")
                    restraint.g1=restraint.g1.difference(atom)
                    logger.info(f"New: {restraint.g1.names}")
                

            
        # check for duplicates in the ligand group g1

        all_restraint_atoms=restraints[0].g1
        for restraint in restraints[1:-1]:
            all_restraint_atoms+=restraint.g1
        logger.info(f"All restraint atoms: {[atom.name for atom in all_restraint_atoms]}")
        duplicate_restraint_atoms=set([atom for atom in all_restraint_atoms if all_restraint_atoms.ix.tolist().count(atom.ix)>1]) # uniquify via set
        
        logger.info(f"Duplicate restraint atoms: {duplicate_restraint_atoms}")

        for atom in duplicate_restraint_atoms:
            assign_atom(atom)
        

    if "manual" in restraint_command_string:
        logger.debug("generating manual selections")
        manual_restraint_list = configuration["simulation"]["manualrestraints"].keys()
        logger.debug(f"Manual restraints defined: {manual_restraint_list}")
        for key in manual_restraint_list:

            restraint = configuration["simulation"]["manualrestraints"][key]
            restraint_kw = {}
            for key in restraint.keys():
                restraint_kw[key] = restraint[key]
            logger.debug(f"Keywords for {restraint}: {restraint_kw}")
            restraints.append(
                Restraint(
                    restraint["group1"], restraint["group2"], universe, **restraint_kw
                )
            )

    return restraints


def write_restraints_yaml(path, system, config, current_step):
    """Takes the full config as read in in utils.py, the information for the intstate and writes the restraints.yaml.

    Args:
        path: the path to write to (e.g. ./combinedstructure/structure/intst2/restraints.yaml
        system: the system object from state.py

    """
    from transformato.mutate import cc_names_struc1, cc_names_struc2

    logger.debug(cc_names_struc1)
    restraints_dict = {
        "intst": {},
        "system": {"structure": {"tlc": f"{system.tlc}"}},
        "simulation": {"restraints": f"{config['simulation']['restraints']}"},
    }

    if "scaling" in config["simulation"]["restraints"]:

        if current_step == 1:
            lambda_value_scaling = 0
        elif current_step == 2:
            lambda_value_scaling = 0.25
        elif current_step == 3:
            lambda_value_scaling = 0.5
        elif current_step == 4:
            lambda_value_scaling = 0.75
        else:
            lambda_value_scaling = 1
    else:
        lambda_value_scaling = 1

    restraints_dict["intst"]["scaling"] = lambda_value_scaling
    if "manual" in config["simulation"]["restraints"]:
        restraints_dict["simulation"]["manualrestraints"] = config["simulation"][
            "manualrestraints"
        ]
    if system.structure == "structure1":
        restraints_dict["system"]["structure"]["ccs"] = cc_names_struc1
    elif system.structure == "structure2":
        restraints_dict["system"]["structure"]["ccs"] = cc_names_struc2

    output = yaml.dump(restraints_dict)
    with open(path, "w") as resyaml:
        resyaml.write(output)
        resyaml.close()
