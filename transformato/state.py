import logging
import os
import shutil

import parmed as pm
import transformato

from .utils import get_toppar_dir
from .mutate import Mutation
from transformato.charmm_factory import charmm_factory, build_reduced_toppar

from typing import List
import types

logger = logging.getLogger(__name__)


class IntermediateStateFactory(object):
    def __init__(self, system: transformato.system, configuration: dict):
        """
        Generate the intermediate directories with for the provided systems with the provided mutations.
        Parameters
        ----------
        system : transformato.system
            definition of the two states for a given system (either waterbox and vacuum or waterbox and complex)
        mutation_list : list
            list of mutations defined by the transformato.ProposeMutationRoute object
        configuration : dict
            configuration dictionary
        """

        self.system = system
        self.path = f"{configuration['system_dir']}/{self.system.name}"
        self._init_base_dir()
        self.configuration = configuration

    def write_state(
        self,
        mutation_conf: List,
        intst_nr: int,
        lambda_value_electrostatic: float = 1.0,
        lambda_value_vdw: float = 1.0,
        common_core_transformation: float = 1.0,
    ):
        """
        write_state Defines everything that is written to the intermediate state directories

        Parameters
        ----------
        mutation_conf : List
            [description]
        intst_nr : int
            [description]
        lambda_value_electrostatic : float, optional
            [description], by default 1.0
        lambda_value_vdw : float, optional
            [description], by default 1.0
        common_core_transformation : float, optional
            [description], by default 1.0

        Returns
        -------
        [type]
            [description]
        """

        logger.info("#########################################")
        logger.info("#########################################")
        output_file_base = self._init_intermediate_state_dir(intst_nr)
        print(f"Writing to {output_file_base}")

        for env in self.system.envs:
            for mutation_type in mutation_conf:

                if (
                    common_core_transformation < 1.0
                ):  # NOTE: THis is inconsisten -- the mutatino_type is the actual mutation in this case
                    mutation_type.mutate(
                        psf=self.system.psfs[env],
                        lambda_value=common_core_transformation,
                    )

                else:
                    mutation_type.print_details()
                    logger.debug(f"Lambda electrostatics: {lambda_value_electrostatic}")
                    logger.debug(f"Lambda vdw: {lambda_value_vdw}")

                    mutator = Mutation(
                        atoms_to_be_mutated=mutation_type.atoms_to_be_mutated,
                        common_core=mutation_type.common_core,
                        dummy_region=mutation_type.dummy_region,
                    )

                    mutator.mutate(
                        psf=self.system.psfs[env],
                        lambda_value_electrostatic=lambda_value_electrostatic,
                        lambda_value_vdw=lambda_value_vdw,
                        vdw_atom_idx=mutation_type.vdw_atom_idx,
                        steric_mutation_to_default=mutation_type.steric_mutation_to_default,
                    )
            self._write_psf(self.system.psfs[env], output_file_base, env)
        self._write_rtf_file(self.system.psfs[env], output_file_base, self.system.tlc)
        self._write_prm_file(self.system.psfs[env], output_file_base, self.system.tlc)
        self._write_toppar_str(
            output_file_base, self.system.tlc
        )  # NOTE for BB: this needs to be adopted for CHARMM
        self._copy_files(output_file_base)  # NOTE FOR BB: this is the main method
        return output_file_base

    def _add_serializer(self, file):
        # adding serializer functions
        f = open(file, "a")
        f.write(
            """
# mw: adding xml serializer to the simulation script
file_name = str(args.psffile).replace('.psf', '')
print(file_name)
serialized_integrator = XmlSerializer.serialize(integrator)
outfile = open(file_name + '_integrator.xml','w')
outfile.write(serialized_integrator)
outfile.close()
serialized_system = XmlSerializer.serialize(system)
outfile = open(file_name + '_system.xml','w')
outfile.write(serialized_system)
outfile.close()
"""
        )
        f.close()

    def _get_simulations_parameters(self):
        prms = {}
        for key in self.configuration["simulation"]["parameters"]:
            prms[key] = self.configuration["simulation"]["parameters"][key]
        return prms

    def _copy_charmm_files(self, intermediate_state_file_path: str):
        """
        _copy_charmm_files Copy CHARMM specific files in running directories

        Parameters
        ----------
        intermediate_state_file_path : [type]
            [description]
        """
        basedir = self.system.charmm_gui_base

        if (
            self.configuration["simulation"]["free-energy-type"]
            == "solvation-free-energy"
        ):

            for env in self.system.envs:
                if env == "waterbox":
                    charmm_input = charmm_factory(
                        configuration=self.configuration,
                        structure=self.system.structure,
                        env="waterbox",
                    )
                    try:
                        with open(
                            f"{intermediate_state_file_path}/charmm_run_waterbox.inp",
                            "w",
                        ) as f:
                            f.write(charmm_input)
                            f.close()
                    except IOError:
                        logger.info(
                            f"Data class could not be created: charmm_lig_in_{env}.inp"
                        )
                    pass
                else:  # vacuum
                    charmm_input = charmm_factory(
                        configuration=self.configuration,
                        structure=self.system.structure,
                        env="vacuum",
                    )
                    try:
                        with open(
                            f"{intermediate_state_file_path}/charmm_run_vacuum.inp", "w"
                        ) as f:
                            f.write(charmm_input)
                            f.close()
                    except IOError:
                        logger.info(
                            f"Data class could not be created: charmm_lig_in_{env}.inp"
                        )
                    pass
                    # CHARMM # NOTE for BB: A CHARMM vacuum script needs to be written -- FOR NOW THIS IS THE OPNEMM script
        elif (
            self.configuration["simulation"]["free-energy-type"]
            == "binding-free-energy"
        ):
            pass
        else:
            pass

        # copy diverse set of helper files for CHARMM
        FILES = [
            "crystal_image.str",
            "step3_pbcsetup.str",
        ]
        for f in FILES:
            try:
                charmm_source = f"{basedir}/waterbox/{f}"
                charmm_target = f"{intermediate_state_file_path}/charmm_{f}"
                shutil.copyfile(charmm_source, charmm_target)
            except OSError:
                logger.critical(f"Could not find file: {f}")
                raise

        # copy crd files
        for env in self.system.envs:
            crd_file_source = f"{basedir}/{env}/{self.configuration['system'][self.system.structure][env]['crd_file_name']}.crd"
            crd_file_target = f"{intermediate_state_file_path}/lig_in_{env}.crd"
            try:
                shutil.copyfile(crd_file_source, crd_file_target)
            except FileNotFoundError:
                logger.warning(
                    f"No crd file found for {env} -- using parmed system structure to create crd file."
                )
                crd_file_target = f"{intermediate_state_file_path}/lig_in_{env}.crd"
                pm.charmm.CharmmCrdFile.write(self.system.psfs[env], crd_file_target)

        # copy rst files
        for env in self.system.envs:
            rst_file_source = f"{basedir}/{env}/{self.configuration['system'][self.system.structure][env]['rst_file_name']}.rst"
            rst_file_target = f"{intermediate_state_file_path}/lig_in_{env}.rst"
            try:
                shutil.copyfile(rst_file_source, rst_file_target)
            except FileNotFoundError:
                logger.warning(
                    f"No restart file found for {env} -- starting simulation from crd file."
                )

        # # copy charmm simulation script
        # charmm_simulation_script_source = f"{basedir}/waterbox/step5_production.inp"
        # charmm_simulation_script_target = (
        #     f"{intermediate_state_file_path}/charmm_production.inp"
        # )
        # shutil.copyfile(
        #     charmm_simulation_script_source, charmm_simulation_script_target
        # )

    def _copy_omm_files(self, intermediate_state_file_path: str):
        """
        _copy_omm_files Copyies the files needed for the production runs with openMM in the intst* directories

        Parameters
        ----------
        intermediate_state_file_path : str
            [description]

        Raises
        ------
        RuntimeError
            [description]
        """
        basedir = self.system.charmm_gui_base

        if (
            self.configuration["simulation"]["free-energy-type"]
            == "solvation-free-energy"
        ):
            # parse omm simulation paramter
            for env in self.system.envs:
                if env == "waterbox":
                    omm_simulation_parameter_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['simulation_parameter']}"
                    omm_simulation_parameter_target = f"{intermediate_state_file_path}/{self.configuration['system'][self.system.structure][env]['intermediate-filename']}"
                    self._overwrite_simulation_script_parameters(
                        omm_simulation_parameter_source, omm_simulation_parameter_target
                    )
                else:  # vacuum
                    used_env = "waterbox"
                    omm_simulation_parameter_source = f"{basedir}/{used_env}/openmm/{self.configuration['system'][self.system.structure][used_env]['simulation_parameter']}"
                    omm_simulation_parameter_target = f"{intermediate_state_file_path}/{self.configuration['system'][self.system.structure][env]['intermediate-filename']}"
                    self._overwrite_simulation_script_parameters(
                        omm_simulation_parameter_source, omm_simulation_parameter_target
                    )

            # copy simulation bash script
            omm_simulation_submit_script_source = (
                f"{self.configuration['bin_dir']}/simulation-solvation-free-energy.sh"
            )
            omm_simulation_submit_script_target = (
                f"{intermediate_state_file_path}/simulation.sh"
            )
            shutil.copyfile(
                omm_simulation_submit_script_source, omm_simulation_submit_script_target
            )
        elif (
            self.configuration["simulation"]["free-energy-type"]
            == "binding-free-energy"
        ):
            # parse omm simulation paramter
            for env in self.system.envs:
                omm_simulation_parameter_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['simulation_parameter']}"
                omm_simulation_parameter_target = f"{intermediate_state_file_path}/{self.configuration['system'][self.system.structure][env]['intermediate-filename']}"
                self._overwrite_simulation_script_parameters(
                    omm_simulation_parameter_source, omm_simulation_parameter_target
                )

            # copy simulation bash script
            omm_simulation_submit_script_source = (
                f"{self.configuration['bin_dir']}/simulation-binding-free-energy.sh"
            )
            omm_simulation_submit_script_target = (
                f"{intermediate_state_file_path}/simulation.sh"
            )
            shutil.copyfile(
                omm_simulation_submit_script_source, omm_simulation_submit_script_target
            )
        else:
            raise RuntimeError(f"Only solvation/binding free energies implemented")

        # copy crd files
        for env in self.system.envs:
            crd_file_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['crd_file_name']}.crd"
            crd_file_target = f"{intermediate_state_file_path}/lig_in_{env}.crd"
            try:
                shutil.copyfile(crd_file_source, crd_file_target)
            except FileNotFoundError:
                logger.warning(
                    f"No crd file found for {env} -- using parmed system structure to create crd file."
                )
                crd_file_target = f"{intermediate_state_file_path}/lig_in_{env}.crd"
                pm.charmm.CharmmCrdFile.write(self.system.psfs[env], crd_file_target)

        # copy rst files
        for env in self.system.envs:
            rst_file_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['rst_file_name']}.rst"
            rst_file_target = f"{intermediate_state_file_path}/lig_in_{env}.rst"
            try:
                shutil.copyfile(rst_file_source, rst_file_target)
            except FileNotFoundError:
                logger.warning(
                    f"No restart file found for {env} -- starting simulation from crd file."
                )

        # copy diverse set of helper functions for openMM
        FILES = [
            "omm_barostat.py",
            "omm_readinputs.py",
            "omm_readparams.py",
            "omm_restraints.py",
            "omm_rewrap.py",
            "omm_vfswitch.py",
            "omm_hmr.py",
        ]
        for f in FILES:
            try:
                omm_source = f"{basedir}/waterbox/openmm/{f}"
                omm_target = f"{intermediate_state_file_path}/{f}"
                shutil.copyfile(omm_source, omm_target)
            except OSError:
                logger.critical(f"Could not find file: {f}")
                raise

        # copy omm simulation script
        omm_simulation_script_source = f"{basedir}/waterbox/openmm/openmm_run.py"
        omm_simulation_script_target = f"{intermediate_state_file_path}/openmm_run.py"
        shutil.copyfile(omm_simulation_script_source, omm_simulation_script_target)
        # add serialization
        self._add_serializer(omm_simulation_script_target)
        omm_simulation_script_source = (
            f"{self.configuration['bin_dir']}/openmm_run_vacuum.py"
        )
        omm_simulation_script_target = (
            f"{intermediate_state_file_path}/openmm_run_vacuum.py"
        )
        shutil.copyfile(omm_simulation_script_source, omm_simulation_script_target)
        self._add_serializer(omm_simulation_script_target)

    def _copy_ligand_specific_top_and_par(
        self, basedir: str, intermediate_state_file_path: str
    ):

        # copy ligand rtf file
        ligand_rtf = f"{basedir}/waterbox/{self.system.tlc.lower()}/{self.system.tlc.lower()}_g.rtf"
        toppar_target = (
            f"{intermediate_state_file_path}/{self.system.tlc.lower()}_g.rtf"
        )
        shutil.copyfile(ligand_rtf, toppar_target)

        # copy ligand prm file
        ligand_prm = f"{basedir}/waterbox/{self.system.tlc.lower()}/{self.system.tlc.lower()}.prm"
        toppar_target = f"{intermediate_state_file_path}/{self.system.tlc.lower()}.prm"
        shutil.copyfile(ligand_prm, toppar_target)

    def _copy_files(
        self, intermediate_state_file_path: str
    ):  # NOTE FOR BB: this is the main thing
        """
        Copy the files from the original CHARMM-GUI output folder in the intermediate directories.
        """

        basedir = self.system.charmm_gui_base

        self._copy_ligand_specific_top_and_par(basedir, intermediate_state_file_path)

        # copy central toppar folder
        toppar_dir = get_toppar_dir()
        toppar_source = f"{toppar_dir}"
        toppar_target = f"{intermediate_state_file_path}/toppar"
        shutil.copytree(toppar_source, toppar_target)

        self._copy_omm_files(intermediate_state_file_path)
        self._copy_charmm_files(intermediate_state_file_path)

    def _overwrite_simulation_script_parameters(
        self, omm_simulation_parameter_source: str, omm_simulation_parameter_target: str
    ):
        """
        _overwrite_simulation_script_parameters Overwrites simulatioo parameters that are defined in omm_simulation_parameter_source

        Parameters
        ----------
        omm_simulation_parameter_source : str
            key:value pair file that defines parameters to overwrite
        omm_simulation_parameter_target : str
            new parameter file for simulation
        """
        overwrite_parameters = self._get_simulations_parameters()

        input_simulation_parameter = open(omm_simulation_parameter_source, "r")
        output_simulation_parameter = open(
            omm_simulation_parameter_target + ".inp", "w+"
        )

        for l in input_simulation_parameter.readlines():
            if l.strip():
                t1, t2 = l.split("=")
                t1 = t1.strip()
                t2, comment = t2.split("#")
                t2 = t2.strip()
                comment = comment.strip()
                if t1 in overwrite_parameters.keys():
                    t2 = overwrite_parameters[t1]
                output_simulation_parameter.write(
                    f"{t1:<25} = {t2:<25} # {comment:<30}\n"
                )
            else:
                output_simulation_parameter.write("\n")
        input_simulation_parameter.close()
        output_simulation_parameter.close()

    def _write_rtf_file(
        self, psf, output_file_base, tlc
    ):  # NOTE: this needs some refactoring!
        """
        Generates the dummy atom parameter rtf.
        """

        header_rtf = """* Dummy atom parameters 
* generated by transformato
*
36  1
"""
        rtf_file_handler = open(output_file_base + "/dummy_atom_definitions.rtf", "w")
        rtf_file_handler.write(header_rtf)
        for atom in psf.view[f":{tlc}"].atoms:
            if hasattr(atom, "initial_type"):
                logger.debug("- Setting dummy parameters ...")
                logger.debug(f"  + Atom-Name: {atom.name}")
                logger.debug(f"  + Atom-Type: {atom.initial_type}")
                logger.debug(f"  + Atom Dummy Type: {atom.type}")

                rtf_file_handler.write(
                    "{:7} {:6} {:6} {:6}\n".format("MASS", "-1", atom.type, atom.mass)
                )
                # TODO: check the  rtf_file_handler.write vs

        rtf_file_handler.close()

    def _write_prm_file(self, psf, output_file_base, tlc):

        header_prm = """* Parameters generated by analogy by
* CHARMM General Force Field (CGenFF) program version 1.0.0
*
! Automatically obtained dummy parameters 
! from transformato
"""

        prm_file_handler = open(f"{output_file_base}/dummy_parameters.prm", "w")
        prm_file_handler.write(header_prm)
        prm_file_handler.write("\nATOMS\n")
        already_seen = list()
        view = psf.view[f":{tlc}"]
        # writing atom parameters
        for atom in view.atoms:
            if hasattr(atom, "initial_type"):
                if set([atom.type]) in already_seen:
                    continue
                else:
                    already_seen.append(set([atom.type]))

                logger.debug("- Setting dummy parameters ...")
                logger.debug(f"  + Atom-Name: {atom.name}")
                logger.debug(f"  + Atom-Type: {atom.initial_type}")
                logger.debug(f"  + Atom Dummy Type: {atom.type}")
                prm_file_handler.write(
                    "{:7} {:6} {:6} {:9.5f}\n".format(
                        "MASS", "-1", atom.type, atom.mass
                    )
                )

        prm_file_handler.write("\n\n")

        ##############################################################################
        # write bond parameters - again there are two ways to use this:
        # - keeping bonded terms between real/dummy and dummy atoms intact
        # - changing bonded parameters between real atoms - this again needs dummy atoms

        prm_file_handler.write("BONDS\n")
        already_seen = []
        for bond in view.bonds:
            atom1, atom2 = bond.atom1, bond.atom2
            if any(hasattr(atom, "initial_type") for atom in [atom1, atom2]):
                if set([atom1.type, atom2.type]) in already_seen:
                    continue
                else:
                    already_seen.append(set([atom1.type, atom2.type]))

                logger.debug(
                    " >> Setting dummy bond parameters for: {} - {}".format(
                        str(atom1.type), str(atom2.type)
                    )
                )
                try:
                    logger.debug(
                        "{:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type),
                            str(atom2.type),
                            bond.mod_type.k,
                            bond.mod_type.req,
                        )
                    )
                    prm_file_handler.write(
                        "{:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type),
                            str(atom2.type),
                            bond.mod_type.k,
                            bond.mod_type.req,
                        )
                    )
                except AttributeError:
                    logger.debug(
                        "{:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type), str(atom2.type), bond.type.k, bond.type.req
                        )
                    )
                    prm_file_handler.write(
                        "{:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type), str(atom2.type), bond.type.k, bond.type.req
                        )
                    )

        #################################################################
        prm_file_handler.write("\n\n")
        prm_file_handler.write("ANGLES\n")
        already_seen = []
        for angle in view.angles:
            atom1, atom2, atom3 = angle.atom1, angle.atom2, angle.atom3
            if any(hasattr(atom, "initial_type") for atom in [atom1, atom2, atom3]):

                if [atom1.type, atom2.type, atom3.type] in already_seen:
                    logger.info(f"Skipping {[atom1.type, atom2.type, atom3.type]}")
                    continue
                else:
                    already_seen.append([atom1.type, atom2.type, atom3.type])

                logger.debug("############################################")
                logger.debug("Printing angle atoms which at least one dummy atom.")
                logger.debug(f"{angle.atom1}, {angle.atom2}, {angle.atom3}")
                logger.debug(
                    f" >> Setting dummy angle parameters for: {atom1.type}-{atom2.type}-{atom3.type}"
                )
                try:
                    prm_file_handler.write(
                        "{:7} {:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type),
                            str(atom2.type),
                            str(atom3.type),
                            angle.mod_type.k,
                            angle.mod_type.theteq,
                        )
                    )
                    logger.debug(
                        "{:7} {:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type),
                            str(atom2.type),
                            str(atom3.type),
                            angle.mod_type.k,
                            angle.mod_type.theteq,
                        )
                    )
                except AttributeError:
                    prm_file_handler.write(
                        "{:7} {:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type),
                            str(atom2.type),
                            str(atom3.type),
                            angle.type.k,
                            angle.type.theteq,
                        )
                    )
                    logger.debug(
                        "{:7} {:7} {:7} {:9.5f} {:9.5f} \n".format(
                            str(atom1.type),
                            str(atom2.type),
                            str(atom3.type),
                            angle.type.k,
                            angle.type.theteq,
                        )
                    )

        #################################################################
        prm_file_handler.write("\n\n")
        prm_file_handler.write("DIHEDRALS\n")
        already_seen = []
        for dihedral in view.dihedrals:
            atom1, atom2, atom3, atom4 = (
                dihedral.atom1,
                dihedral.atom2,
                dihedral.atom3,
                dihedral.atom4,
            )
            if any(
                hasattr(atom, "initial_type") for atom in [atom1, atom2, atom3, atom4]
            ):
                if [atom1.type, atom2.type, atom3.type, atom4.type] in already_seen:
                    continue
                else:
                    already_seen.append(
                        [atom1.type, atom2.type, atom3.type, atom4.type]
                    )

                logger.debug(
                    f" >> Setting dummy dihedral parameters for: {atom1.type}-{atom2.type}-{atom3.type}-{atom4.type}"
                )
                try:
                    for dihedral_type in dihedral.mod_type:
                        prm_file_handler.write(
                            "{:7} {:7} {:7} {:7} {:6.5f} {:9.5f} {:9.5f} \n".format(
                                str(atom1.type),
                                str(atom2.type),
                                str(atom3.type),
                                str(atom4.type),
                                dihedral_type.phi_k,
                                dihedral_type.per,
                                dihedral_type.phase,
                            )
                        )
                except AttributeError:
                    for dihedral_type in dihedral.type:
                        prm_file_handler.write(
                            "{:7} {:7} {:7} {:7} {:6.5f} {:9.5f} {:9.5f} \n".format(
                                str(atom1.type),
                                str(atom2.type),
                                str(atom3.type),
                                str(atom4.type),
                                dihedral_type.phi_k,
                                dihedral_type.per,
                                dihedral_type.phase,
                            )
                        )

        #################################################################
        # get all unique improper and parameters
        prm_file_handler.write("\n\n")
        prm_file_handler.write("IMPROPERS\n")
        for impr in view.impropers:
            atom1, atom2, atom3, atom4 = impr.atom1, impr.atom2, impr.atom3, impr.atom4
            if any(
                hasattr(atom, "initial_type") for atom in [atom1, atom2, atom3, atom4]
            ):
                # print('>> Setting dummy improper parameters for: {}-{}-{}-{}'.format(str(atom1.type),str(atom2.type),str(atom3.type),str(atom4.type)))
                # carefull with this solution - > central atom has to be set in the beginning
                prm_file_handler.write(
                    "{:7} {:7} {:7} {:7} {:9.5f} {:9.5f} \n".format(
                        str(atom1.type),
                        str(atom2.type),
                        str(atom3.type),
                        str(atom4.type),
                        impr.type.psi_k,
                        impr.type.psi_eq,
                    )
                )

        #################################################################
        prm_file_handler.write("\n\n")
        prm_file_handler.write(
            """NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5"""
        )
        prm_file_handler.write("\n\n")

        for atom in view.atoms:
            if hasattr(atom, "initial_type"):
                try:
                    prm_file_handler.write(
                        "{:7} {:6} {:9.5f} {:9.5f}\n".format(
                            atom.type, 0.0, atom.mod_type.epsilon, atom.mod_type.rmin
                        )
                    )
                except AttributeError:
                    prm_file_handler.write(
                        "{:7} {:6} {:9.5f} {:9.5f}\n".format(
                            atom.type, 0.0, atom.epsilon, atom.rmin
                        )
                    )

        prm_file_handler.write("\n")
        prm_file_handler.write("END")
        prm_file_handler.close()

    def _init_base_dir(self):
        """
        Generates the base directory which all intermediate states are located.
        """

        if os.path.isdir(self.path):
            shutil.rmtree(self.path)
            os.makedirs(self.path)
        else:
            os.makedirs(self.path)

    def _write_toppar_str(self, output_file_base, tlc):

        toppar_format = f"""
toppar/top_all36_prot.rtf
toppar/par_all36m_prot.prm
toppar/top_all36_na.rtf
toppar/par_all36_na.prm
toppar/top_all36_carb.rtf
toppar/par_all36_carb.prm
toppar/top_all36_lipid.rtf
toppar/par_all36_lipid.prm
toppar/top_all36_cgenff.rtf
toppar/par_all36_cgenff.prm
toppar/toppar_water_ions.str
toppar/toppar_dum_noble_gases.str
toppar/toppar_all36_prot_c36m_d_aminoacids.str
toppar/toppar_all36_prot_fluoro_alkanes.str
toppar/toppar_all36_prot_heme.str
toppar/toppar_all36_prot_na_combined.str
toppar/toppar_all36_prot_retinol.str
toppar/toppar_all36_na_nad_ppi.str
toppar/toppar_all36_lipid_bacterial.str
toppar/toppar_all36_lipid_cardiolipin.str
toppar/toppar_all36_lipid_cholesterol.str
toppar/toppar_all36_lipid_inositol.str
toppar/toppar_all36_lipid_lps.str
toppar/toppar_all36_lipid_miscellaneous.str
toppar/toppar_all36_lipid_model.str
toppar/toppar_all36_lipid_prot.str
toppar/toppar_all36_lipid_sphingo.str
{tlc.lower()}_g.rtf
{tlc.lower()}.prm
dummy_atom_definitions.rtf
dummy_parameters.prm
"""

        f = open(f"{output_file_base}/toppar.str", "w+")
        f.write(toppar_format)
        f.close()

        # TODO for BB: write charmm toppar.str file: charmm_toppar.str
        charmm_toppar = build_reduced_toppar(tlc.lower())

        f = open(f"{output_file_base}/charmm_toppar.str", "w+")
        f.write(charmm_toppar)
        f.close()

    def _write_psf(self, psf, output_file_base: str, env: str):
        """
        Writes the new psf.
        """

        psf.write_psf(f"{output_file_base}/lig_in_{env}.psf")
        psf.write_pdb(f"{output_file_base}/lig_in_{env}.pdb")

    def _init_intermediate_state_dir(self, nr: int):
        """
        Generates the intermediate state directory.
        """
        output_file_base = f"{self.path}/intst{nr}/"

        logger.info(f" - Created directory: - {os.path.abspath(output_file_base)}")
        os.makedirs(output_file_base)
        logger.info(f" - Writing in - {os.path.abspath(output_file_base)}")
        return output_file_base
