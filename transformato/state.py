import logging
import os
import shutil
import glob
from io import StringIO
from typing import List

import parmed as pm

import transformato
from transformato.charmm_factory import CharmmFactory
from transformato.mutate import Mutation
from transformato.restraints import write_restraints_yaml
from transformato.utils import get_toppar_dir, psf_correction, check_switching_function

logger = logging.getLogger(__name__)


class IntermediateStateFactory(object):
    def __init__(
        self,
        system: transformato.system.SystemStructure,
        configuration: dict,
        multiple_runs=False,
    ):
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
        multiple_runs : bool or int
            defines how many runs should be performed per intst state, if not used it is 0 (=False) else give an integer number
        """

        self.system = system
        self.path = f"{configuration['system_dir']}/{self.system.name}"
        self.configuration = configuration
        self._init_base_dir()
        self.vdw_switch: str
        self.charmm_factory = CharmmFactory(configuration, self.system.structure)
        self.output_files = []
        self.current_step = 1
        self.multiple_runs = multiple_runs

    def endstate_correction(self):

        logger.info(f"Will create script for endstate correction")
        try:
            os.makedirs(f"{self.path}/../endstate_correction")
        except:
            logger.info(f"Folder for endstate correction exist already")

        # copy submit script to the newly created folder
        submit_switching_script_source = (
            f"{self.configuration['bin_dir']}/slurm_switching.sh"
        )
        submit_switchting_script_target = (
            f"{self.path}/../endstate_correction/slurm_switching.sh"
        )
        shutil.copyfile(submit_switching_script_source, submit_switchting_script_target)

        # modify the perform_correction file from transformato bin and save it in the endstate_correcition folder
        endstate_correction_script_source = (
            f"{self.configuration['bin_dir']}/perform_correction.py"
        )
        endstate_correction_script_target = (
            f"{self.path}/../endstate_correction/perform_correction.py"
        )
        fin = open(endstate_correction_script_source, "rt")
        fout = open(endstate_correction_script_target, "wt")

        for line in fin:
            if "NAMEofSYSTEM" in line:
                fout.write(
                    line.replace(
                        "NAMEofSYSTEM",
                        f'"{self.configuration["system"]["structure1"]["name"]}"',
                    )
                )
            elif "TLC" in line:
                fout.write(
                    line.replace(
                        "TLC", f'"{self.configuration["system"]["structure1"]["tlc"]}"'
                    )
                )
            else:
                fout.write(line)
        fin.close()
        fout.close()

    def write_state(
        self,
        mutation_conf: List,
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
        lambda_value_electrostatic : float, optional
            [description], by default 1.0
        lambda_value_vdw : float, optional
            [description], by default 1.0
        common_core_transformation : float, optional
            [description], by default 1.0

        """

        logger.debug("#########################################")
        output_file_base = self._init_intermediate_state_dir(self.current_step)
        logger.info(f"Writing to {output_file_base}")

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

                    logger.debug(f"Lambda electrostatics: {lambda_value_electrostatic}")
                    logger.debug(f"Lambda vdw: {lambda_value_vdw}")

                    mutator = Mutation(
                        atoms_to_be_mutated=mutation_type.atoms_to_be_mutated,
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
        self._write_toppar_str(output_file_base)
        self._copy_files(output_file_base)

        # Create run folder for dcd output for each intst state
        if self.multiple_runs:
            if type(self.multiple_runs) == str:
                os.makedirs(f"{output_file_base}/{self.multiple_runs}")
            else:
                for i in range(1, self.multiple_runs + 1):
                    os.makedirs(f"{output_file_base}/run_{i}")

        self.output_files.append(output_file_base)

        # Used for restraints:
        if "restraints" in self.configuration["simulation"].keys():

            logger.info(
                "Found restraints in configuration file - writing restraints.yaml"
            )
            write_restraints_yaml(
                f"{output_file_base}/restraints.yaml",
                self.system,
                self.configuration,
                self.current_step,
            )

        self.current_step += 1

    def _add_serializer(self, file):
        # adding serializer functions
        with open(file, "a") as f:
            f.write(
                """
# mw: adding xml serializer to the simulation script
file_name = str(args.psffile).replace('.psf', '')
print(file_name)
with open(file_name + '_integrator.xml','w') as outfile:
    outfile.write(XmlSerializer.serialize(integrator))
with open(file_name + '_system.xml','w') as outfile:
    outfile.write(XmlSerializer.serialize(system))
"""
            )

    def _get_simulations_parameters(self):
        prms = {}
        for key in self.configuration["simulation"]["parameters"]:
            prms[key] = self.configuration["simulation"]["parameters"][key]
        return prms

    def _write_workload_preamble(self, filepath):
        """
        Prepends the preamble of the selected workload manager in the file specified by filepath

        The preamble is the *-preamble.sh file found in transformato/bin

        Args:
            filepath (str): Path to the file (usually _script_target)

        """

        workloadmanager = self.configuration["simulation"]["workload-manager"]
        with open(
            f"{self.configuration['bin_dir']}/{workloadmanager}-preamble.sh", "r"
        ) as preamblefile:
            preamble = preamblefile.read()

        with open(filepath, "r+") as script_to_prepend_to:
            content = script_to_prepend_to.read()
            script_to_prepend_to.seek(0)
            script_to_prepend_to.write(f"{preamble}\n{content}")

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
            self.configuration["simulation"]["free-energy-type"] == "rsfe"
            or self.configuration["simulation"]["free-energy-type"] == "asfe"
        ):
            # copy simulation bash script
            charmm_simulation_submit_script_source = (
                f"{self.configuration['bin_dir']}/simulation-rsfe_charmm.sh"
            )
            charmm_simulation_submit_script_target = (
                f"{intermediate_state_file_path}/simulation_charmm.sh"
            )
            shutil.copyfile(
                charmm_simulation_submit_script_source,
                charmm_simulation_submit_script_target,
            )

            for env in self.system.envs:
                if env == "waterbox":
                    # write charmm production scripte
                    charmm_input = self.charmm_factory.generate_CHARMM_production_files(
                        env,
                    )
                    with open(
                        f"{intermediate_state_file_path}/charmm_run_waterbox.inp", "w+"
                    ) as f:
                        f.write(charmm_input)

                    # write charmm postproduction script
                    charmm_input = (
                        self.charmm_factory.generate_CHARMM_postprocessing_files(
                            env,
                        )
                    )
                    with open(
                        f"{intermediate_state_file_path}/charmm_evaluate_energy_in_{env}.inp",
                        "w+",
                    ) as f:
                        f.write(charmm_input)

                elif env == "vacuum":  # vacuum
                    # write charmm production scripte
                    charmm_input = self.charmm_factory.generate_CHARMM_production_files(
                        env,
                    )
                    with open(
                        f"{intermediate_state_file_path}/charmm_run_vacuum.inp", "w"
                    ) as f:
                        f.write(charmm_input)
                    # write charmm postproduction script
                    charmm_input = (
                        self.charmm_factory.generate_CHARMM_postprocessing_files(
                            env,
                        )
                    )
                    with open(
                        f"{intermediate_state_file_path}/charmm_evaluate_energy_in_{env}.inp",
                        "w+",
                    ) as f:
                        f.write(charmm_input)

                else:
                    raise NotImplementedError()

        elif self.configuration["simulation"]["free-energy-type"] == "rbfe":
            # copy simulation bash script
            charmm_simulation_submit_script_source = (
                f"{self.configuration['bin_dir']}/simulation-rbfe_charmm.sh"
            )
            charmm_simulation_submit_script_target = (
                f"{intermediate_state_file_path}/simulation_charmm.sh"
            )
            shutil.copyfile(
                charmm_simulation_submit_script_source,
                charmm_simulation_submit_script_target,
            )

            for env in self.system.envs:
                if env == "waterbox" or env == "complex":
                    # write charmm production scripte
                    charmm_input = self.charmm_factory.generate_CHARMM_production_files(
                        env,
                    )
                    with open(
                        f"{intermediate_state_file_path}/charmm_run_{env}.inp", "w+"
                    ) as f:
                        f.write(charmm_input)

                    # write charmm postproduction script
                    charmm_input = (
                        self.charmm_factory.generate_CHARMM_postprocessing_files(
                            env,
                        )
                    )
                    with open(
                        f"{intermediate_state_file_path}/charmm_evaluate_energy_in_{env}.inp",
                        "w+",
                    ) as f:
                        f.write(charmm_input)

                else:
                    raise NotImplementedError()

        else:
            pass

        # Prepend workload manager instructions
        self._write_workload_preamble(charmm_simulation_submit_script_target)

        # copy diverse set of helper files for CHARMM
        for env in self.system.envs:
            if env != "vacuum":
                FILES = [
                    "crystal_image.str",
                    "step3_pbcsetup.str",
                ]
                for f in FILES:
                    try:
                        charmm_source = f"{basedir}/{env}/{f}"
                        charmm_target = (
                            f"{intermediate_state_file_path}/charmm_{env}_{f}"
                        )
                        shutil.copyfile(charmm_source, charmm_target)
                    except FileNotFoundError:
                        logger.critical(f"Could not find file: {f}")
                        raise

            # copy rst files
            rst_file_source = f"{basedir}/{env}/{self.configuration['system'][self.system.structure][env]['rst_file_name']}.rst"
            rst_file_target = f"{intermediate_state_file_path}/lig_in_{env}.rst"
            try:
                shutil.copyfile(rst_file_source, rst_file_target)
            except FileNotFoundError:
                logger.warning(
                    f"No restart file found for {env} -- starting simulation from crd file."
                )

    def _modify_submit_script(self, shFile):
        # submit script is modified to loop over all runs depending on the number defined in the consecutive
        # runs
        logger.info(f"We will manipulate the submit script for use with multiple runs")

        fout = open(f"{shFile}.tmp", "wt")
        with open(f"{shFile}", "r+") as f:
            for line in f:
                if line.startswith(f"input=lig_in_"):
                    fout.write("for i in {1.." + f"{self.multiple_runs}" + "};\n")
                    fout.write("do \n")
                    fout.write(line)
                elif line.startswith("python -u openmm"):
                    line = line.replace("${istep}.dcd", "run_${i}/${istep}.dcd")
                    fout.write(
                        line.replace(
                            line.split()[-1], "run_${i}/" + f"{line.split()[-1]}"
                        )
                    )
                    fout.write(f"done \n")
                else:
                    fout.write(line)
        fout.close()
        shutil.move(fout.name, shFile)

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
            self.configuration["simulation"]["free-energy-type"] == "rsfe"
            or self.configuration["simulation"]["free-energy-type"] == "asfe"
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
                f"{self.configuration['bin_dir']}/simulation-rsfe.sh"
            )
            omm_simulation_submit_script_target = (
                f"{intermediate_state_file_path}/simulation.sh"
            )
            shutil.copyfile(
                omm_simulation_submit_script_source, omm_simulation_submit_script_target
            )

        elif self.configuration["simulation"]["free-energy-type"] == "rbfe":
            # parse omm simulation paramter
            for env in self.system.envs:
                omm_simulation_parameter_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['simulation_parameter']}"
                omm_simulation_parameter_target = f"{intermediate_state_file_path}/{self.configuration['system'][self.system.structure][env]['intermediate-filename']}"
                self._overwrite_simulation_script_parameters(
                    omm_simulation_parameter_source, omm_simulation_parameter_target
                )

            # copy simulation bash script
            omm_simulation_submit_script_source = (
                f"{self.configuration['bin_dir']}/simulation-rbfe.sh"
            )
            omm_simulation_submit_script_target = (
                f"{intermediate_state_file_path}/simulation.sh"
            )
            shutil.copyfile(
                omm_simulation_submit_script_source, omm_simulation_submit_script_target
            )
        else:
            raise RuntimeError(f"Only solvation/binding free energies implemented")

        if self.multiple_runs:
            self._modify_submit_script(omm_simulation_submit_script_target)

        # Prepend workload manager instructions
        self._write_workload_preamble(omm_simulation_submit_script_target)

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

        # copy omm simulation script
        # start with waterbox
        omm_simulation_script_source = f"{self.configuration['bin_dir']}/openmm_run.py"
        omm_simulation_script_target = f"{intermediate_state_file_path}/openmm_run.py"
        shutil.copyfile(omm_simulation_script_source, omm_simulation_script_target)
        # add serialization
        self._check_hmr(omm_simulation_script_target)
        self._add_serializer(omm_simulation_script_target)
        self._change_platform(omm_simulation_script_target)
        check_switching_function(self.vdw_switch)

        if (
            self.configuration["simulation"]["free-energy-type"] == "rsfe"
            or self.configuration["simulation"]["free-energy-type"] == "asfe"
        ):
            # add vacuum scripts
            omm_simulation_script_source = (
                f"{self.configuration['bin_dir']}/openmm_run_vacuum.py"
            )
            omm_simulation_script_target = (
                f"{intermediate_state_file_path}/openmm_run_vacuum.py"
            )
            shutil.copyfile(omm_simulation_script_source, omm_simulation_script_target)
            self._check_hmr(omm_simulation_script_target)
            self._add_serializer(omm_simulation_script_target)
            self._change_platform(omm_simulation_script_target)

    def _check_hmr(self, file):
        "add hmr if requested in config file"
        if self.configuration["simulation"].get("HMR", False):
            with open(file, "r") as f:
                with open(f"{file}_tmp", "w+") as g:
                    for line in f.readlines():
                        if "system = psf.createSystem" in line:
                            g.write("# Adding HMR\n")
                            g.write(
                                "nboptions['hydrogenMass'] = 1.5 *atom_mass_units\n"
                            )
                            g.write(line)
                        else:
                            g.write(line)
            shutil.move(f"{file}_tmp", file)

    def _change_platform(self, file):
        # changin input script
        f = open(file, "r")
        g = open(f"{file}_tmp", "w+")
        i = 0  # counting lines

        if self.configuration["simulation"]["GPU"]:
            logger.debug("Preparing for CUDA")
            for line in f.readlines():
                if "# Set platform" in line and i == 0:
                    i += 1
                    g.write(line)
                elif i == 1:
                    i += 1
                    g.write("platform = Platform.getPlatformByName('CUDA')\n")
                elif i == 2:
                    i += 2
                    g.write("prop = dict(CudaPrecision='mixed')\n")
                else:
                    g.write(line)
        else:
            for line in f.readlines():
                if "# Set platform" in line and i == 0:
                    i += 1
                    g.write(line)
                elif i == 1:
                    i += 1
                    g.write("platform = Platform.getPlatformByName('CPU')\n")
                elif i == 2:
                    i += 2
                    g.write("prop = dict()\n")
                else:
                    g.write(line)

        f.close()
        g.close()
        shutil.move(f"{file}_tmp", file)

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

    def _copy_ligand_specific_str(
        self, basedir: str, intermediate_state_file_path: str
    ):

        # copy ligand rtf file
        ligand_rtf = f"{basedir}/waterbox/{self.system.tlc.lower()}/{self.system.tlc.lower()}.str"
        toppar_target = f"{intermediate_state_file_path}/{self.system.tlc.lower()}.str"
        shutil.copyfile(ligand_rtf, toppar_target)

    def _copy_crd_file(self, intermediate_state_file_path: str):

        basedir = self.system.charmm_gui_base
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

    def _copy_files(self, intermediate_state_file_path: str):
        """
        Copy the files from the original CHARMM-GUI output folder in the intermediate directories.
        """

        basedir = self.system.charmm_gui_base

        if os.path.exists(f"{basedir}/waterbox/{self.system.tlc.lower()}/"):
            try:
                self._copy_ligand_specific_top_and_par(
                    basedir, intermediate_state_file_path
                )
            except:
                self._copy_ligand_specific_str(basedir, intermediate_state_file_path)

        # copy crd file
        self._copy_crd_file((intermediate_state_file_path))

        # copy openMM and charmm specific scripts
        self._copy_omm_files(intermediate_state_file_path)
        self._copy_charmm_files(intermediate_state_file_path)

    def _overwrite_simulation_script_parameters(
        self, omm_simulation_parameter_source: str, omm_simulation_parameter_target: str
    ):
        """
        _overwrite_simulation_script_parameters changes parameters that are defined in omm_simulation_parameter_source

        Parameters
        ----------
        omm_simulation_parameter_source : str
            key:value pair file that defines parameters to overwrite
        omm_simulation_parameter_target : str
            new parameter file for simulation
        """
        import copy

        overwrite_parameters = copy.deepcopy(self._get_simulations_parameters())

        if "vdw" in overwrite_parameters.keys():
            if overwrite_parameters["vdw"].lower() not in [
                "force-switch",
                "no-switch",
                "switch",
            ]:
                raise NotImplementedError(
                    f"switch: {overwrite_parameters['vdw']} not implemented."
                )
            self.vdw_switch = overwrite_parameters["vdw"].lower()
        else:
            self.vdw_switch = "force-switch"  # default for now
            logger.critical("Setting switching function to vfswitch")

        input_simulation_parameter = open(omm_simulation_parameter_source, "r")
        output_simulation_parameter = open(
            omm_simulation_parameter_target + ".inp", "w+"
        )

        common_keywords = [
            "nstep",
            "nstdcd",
            "nstout",
            "cons",
            "dt",
            "mini_nstep",
        ]

        # test that common keywords are in yaml files
        if not all(elem in overwrite_parameters.keys() for elem in common_keywords):
            for elem in common_keywords:
                if elem not in overwrite_parameters.keys():
                    if elem == "cons":
                        logger.critical(f"###################")
                        logger.critical(
                            f"{elem} is not set in config yaml. This is very likely a mistake!"
                        )
                        logger.critical(f"###################")
                    else:
                        logger.critical(f"###################")
                        logger.critical(
                            f"{elem} is not set in config yaml. Was this a mistake?"
                        )
                        logger.critical(f"###################")

        for l in input_simulation_parameter.readlines():
            if l.strip():
                t1, t2_comment = [e.strip() for e in l.split("=")]
                t2, comment = [e.strip() for e in t2_comment.split("#")]
                comment = comment.strip()
                if t1 in overwrite_parameters.keys():
                    t2 = overwrite_parameters[t1]
                    del overwrite_parameters[t1]  # remove from dict
                if t1 == "vdw":
                    t2 = t2.capitalize()
                output_simulation_parameter.write(
                    f"{t1:<25} = {t2:<25} # {comment:<30}\n"
                )
            else:
                output_simulation_parameter.write("\n")

        # set parameters that have no equivalent in the pregenerated parameter file
        for t1 in overwrite_parameters.keys():
            t2 = overwrite_parameters[t1]
            output_simulation_parameter.write(
                f"{t1:<25} = {t2:<25} # some new options\n"
            )

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
                    logger.debug(f"Skipping {[atom1.type, atom2.type, atom3.type]}")
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
        Generates the base directory in which all intermediate states are located
        and create the central toppar dir.
        """

        if os.path.isdir(self.path):
            shutil.rmtree(self.path)
            os.makedirs(self.path)
        else:
            os.makedirs(self.path)

        # copy central toppar folder

        basedir = self.system.charmm_gui_base
        toppar_source = f"{basedir}/waterbox/toppar"

        if os.path.isdir(toppar_source):
            pass
        else:
            toppar_dir = get_toppar_dir()
            toppar_source = f"{toppar_dir}"
            logger.warning(
                f"will use the toppar version available within the Transformato package"
            )

        toppar_target = f"{self.configuration['system_dir']}/toppar"
        shutil.copytree(toppar_source, toppar_target, dirs_exist_ok=True)

    def _write_toppar_str(self, output_file_base):

        from transformato.system import parameter_files

        base = "../../toppar"
        toppar_format = ""
        # Only necessary because the variable self.configuration["system"]["structure2"]["tlc"] gives a keyError for asfe
        # since its not defined there
        try:
            if self.configuration["system"]["structure2"]["tlc"].lower():
                struct2 = self.configuration["system"]["structure2"]["tlc"].lower()
            else:
                struct2 = "emptyString"

        except KeyError:
            struct2 = "emptyString"

        for i in parameter_files:
            # the parameterfiles for customized ligands are in another directory and are considered below
            if self.configuration["system"]["structure1"]["tlc"].lower() in i:
                pass
            if struct2 in i:
                pass
            else:
                toppar_format += f"{base}/{i.split('/')[-1]} \n"

        if os.path.isfile(
            f"{self.system.charmm_gui_base}/waterbox/{self.system.tlc.lower()}/{self.system.tlc.lower()}_g.rtf"
        ):
            toppar_format += f"""{self.system.tlc.lower()}_g.rtf
{self.system.tlc.lower()}.prm
dummy_atom_definitions.rtf
dummy_parameters.prm
"""
        elif os.path.isfile(
            f"{self.system.charmm_gui_base}/waterbox/{self.system.tlc.lower()}/{self.system.tlc.lower()}.str"
        ):
            toppar_format += f"""{self.system.tlc.lower()}.str
dummy_atom_definitions.rtf
dummy_parameters.prm
"""
        else:
            toppar_format += f"""dummy_atom_definitions.rtf
dummy_parameters.prm
"""

        f = open(f"{output_file_base}/toppar.str", "w+")
        f.write(toppar_format)
        f.close()

        # write charmm_toppar.str
        charmm_toppar = self.charmm_factory.build_reduced_toppar(
            self.system.tlc.lower(), self.system.charmm_gui_base
        )

        f = open(f"{output_file_base}/charmm_toppar.str", "w+")
        f.write(charmm_toppar)
        f.close()

    @staticmethod
    def _write_psf(psf, output_file_base: str, env: str):
        """
        Writes the new psf and pdb file.
        """

        with open(f"{output_file_base}/lig_in_{env}.psf", "w+") as f:
            psf.write_psf(f)

        string_object = StringIO()
        psf.write_psf(string_object)
        # read in psf and correct some aspects of the file not suitable for CHARMM
        corrected_psf = psf_correction(string_object)
        with open(f"{output_file_base}/lig_in_{env}.psf", "w+") as f:
            f.write(corrected_psf)
        # write pdb
        psf.write_pdb(f"{output_file_base}/lig_in_{env}.pdb")

    def _init_intermediate_state_dir(self, nr: int):
        """
        Generates the intermediate state directory.
        """
        output_file_base = f"{self.path}/intst{nr}/"

        logger.debug(f" - Created directory: - {os.path.abspath(output_file_base)}")
        os.makedirs(output_file_base)
        logger.info(f" - Writing in - {os.path.abspath(output_file_base)}")
        return output_file_base
