import logging
import os
import glob
import re
from collections import defaultdict
from typing import Tuple

import networkx as nx
import parmed as pm
from rdkit import Chem

from transformato.utils import get_toppar_dir

logger = logging.getLogger(__name__)


class SystemStructure(object):
    def __init__(self, configuration: dict, structure: str):
        """
        A class that contains all informations for a single small ligand in different environments.
        Everything is constructed automatically from the charmm-gui folder.
        Parameters
        ----------
        configuration: dict
            the configuration dictionary obtained with utils.load_config_yaml
        structure: str
            either 'structure1' or 'structure2'
        """

        self.structure: str = structure
        self.name: str = configuration["system"][structure]["name"]
        self.tlc: str = configuration["system"][structure]["tlc"]
        self.charmm_gui_base: str = configuration["system"][structure]["charmm_gui_dir"]
        self.psfs: defaultdict = defaultdict(pm.charmm.CharmmPsfFile)
        self.offset: defaultdict = defaultdict(int)
        # self.parameter = self._read_parameters("waterbox")
        self.cgenff_version: float
        self.envs = set()
        # running a binding-free energy calculation?
        if configuration["simulation"]["free-energy-type"] == "rbfe":
            self.envs = set(["complex", "waterbox"])
        elif (
            configuration["simulation"]["free-energy-type"] == "rsfe"
            or configuration["simulation"]["free-energy-type"] == "asfe"
        ):
            self.envs = set(["waterbox", "vacuum"])
        else:
            raise NotImplementedError(
                "only binding and solvation free energy implemented."
            )

        for env in self.envs:
            # set up system
            self.psfs[env] = self._initialize_system(configuration, env)
            # load parameters, by default only the most important toppar files are loaded
            try:
                parameter = self._read_parameters(env)
                self.psfs[env].load_parameters(parameter)
            except pm.exceptions.ParameterError:
                parameter = self._read_parameters(env, full_set=True)
                self.psfs[env].load_parameters(parameter)

            # for point mutation, expects TIP3P water and NaCl as ions
            if not self.tlc:
                self.tlc = str(self._get_tlc(self.psfs[env]))

            # get offset
            self.offset[env] = self._determine_offset_and_set_possible_dummy_properties(
                self.psfs[env]
            )

        # generate rdkit mol object of small molecule
        self.mol: Chem.Mol = self._generate_rdkit_mol(
            "waterbox", self.psfs["waterbox"][f":{self.tlc}"]
        )  # for RBFE this was called "complex"
        self.graph: nx.Graph = self.mol_to_nx(self.mol)

        # For RBFE:
        # generate rdkit mol object of small molecule
        # self.mol: Chem.Mol = self._generate_rdkit_mol(
        #     "complex", self.psfs["complex"][f":{self.tlc}"]
        # )
        # self.graph: nx.Graph = self.mol_to_nx(self.mol)

    def _get_tlc(self, psf) -> str:
        """
        If no information about a small ligend is available this function will try to find
        all residues which are in the first chain
        """

        # This works only if TIP3 water and NaCl as ions is used
        # it will consider other ions as residue as well
        # It first checks for the chains and will always take the FIRST
        # chain
        chains = []
        for atom in psf.view:
            if atom.residue.name not in chains:
                chains.append(atom.residue.chain)

        tlc = ""
        for atom in psf.view["!(:TIP3:SOD:POT:CLA)"]:
            if atom.residue.name not in tlc and atom.residue.chain == chains[0]:
                tlc += f":{atom.residue.name}"

        assert len(tlc) < 20

        # Remove first colon so it looks e.g. like this: GUA:CYT:URA
        return tlc[1:]

    def _set_hmr(self, configuration: dict, env: str):
        # check if HMR is set
        if configuration["simulation"].get("HMR", False):
            logger.info("Using HMR.")
            for bond in self.psfs[env].bonds:
                # Only take the ones with at least one hydrogen
                atom1, atom2 = bond.atom1, bond.atom2
                # print(dir(atom1))
                if atom1.atomic_number != 1 and atom2.atomic_number != 1:
                    continue
                if "TIP3" in [atom1.residue.name, atom2.residue.name]:
                    continue
                if atom2.atomic_number == 1:
                    atom1, atom2 = (
                        atom2,
                        atom1,
                    )  # now atom1 is hydrogen for sure
                if atom2.atomic_number != 1:
                    mass1 = atom1.mass
                    mass2 = atom2.mass
                    new_mass1 = mass1 * 3.0
                    new_mass2 = mass2 - new_mass1 + mass1
                    atom1.mass = new_mass1
                    atom2.mass = new_mass2

    @staticmethod
    def mol_to_nx(mol: Chem.Mol):
        try:
            from tf_routes import preprocessing

            return preprocessing._mol_to_nx_full_weight(mol)
        except ModuleNotFoundError:

            G = nx.Graph()

            for atom in mol.GetAtoms():
                G.add_node(
                    atom.GetIdx(),
                    atomic_num=atom.GetAtomicNum(),
                    formal_charge=atom.GetFormalCharge(),
                    chiral_tag=atom.GetChiralTag(),
                    hybridization=atom.GetHybridization(),
                    num_explicit_hs=atom.GetNumExplicitHs(),
                    is_aromatic=atom.GetIsAromatic(),
                )

            for bond in mol.GetBonds():
                G.add_edge(
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    bond_type=bond.GetBondType(),
                )

            return G

    def _read_parameters(
        self, env: str, full_set: bool = False
    ) -> pm.charmm.CharmmParameterSet:
        """
        Reads in topparameters from a toppar dir and ligand specific parameters.
        Parameters
        ----------
        env: str
            waterbox,complex or vacuum
        full_set: bool
            wheather all files mentioned in the openmm/toppar.str should be read in, default is False
        Returns
        ----------
        parameter : pm.charmm.CharmmParameterSet
            parameters obtained from the CHARMM-GUI output dir.
        """
        # save list of files for creation of toppar.str
        global parameter_files
        parameter_files = tuple()

        # the parameters for the vacuum system is parsed from the waterbox charmm-gui directory
        if env == "vacuum":
            env = "waterbox"

        charmm_gui_env = self.charmm_gui_base + env
        tlc_lower = str(self.tlc).lower()
        toppar_dir = f"{charmm_gui_env}/toppar"

        # check if toppar dir is available in CHARMM-GUI folder, if not fall back to toppar dir from transformato
        if os.path.isdir(toppar_dir):
            logger.info(
                f"Using the toppar directory from the CHARMM-GUI folder; {toppar_dir}"
            )
        else:
            toppar_dir = get_toppar_dir()
            logger.info(
                f"Using the toppar directory provided in the transformato package; {toppar_dir}"
            )

        # if custom parameter are added they are located in l1,l2
        l1 = f"{charmm_gui_env}/{tlc_lower}/{tlc_lower}.rtf"
        l2 = f"{charmm_gui_env}/{tlc_lower}/{tlc_lower}.prm"
        l3 = f"{charmm_gui_env}/{tlc_lower}/{tlc_lower}.str"

        for file_path in [l1, l2, l3]:
            if os.path.isfile(file_path):
                parameter_files += (file_path,)
            else:
                logger.critical(
                    f"Custom ligand parameters are not present in {file_path}"
                )

        # check cgenff versions
        if parameter_files:
            with open(parameter_files[0]) as f:
                _ = f.readline()
                cgenff = f.readline().rstrip()
                logger.info(f"CGenFF version: {cgenff}")
                cgenff_version = re.findall("\d+\.\d+", cgenff)[0]
                self.cgenff_version = float(cgenff_version)

        if full_set:
            with open(f"{charmm_gui_env}/openmm/toppar.str", "r") as ommtopparstream:
                for line in ommtopparstream:
                    if tlc_lower:
                        if line.strip() != "" and not tlc_lower in line:
                            filename = line.strip("\n").split("/")[-1]
                            parameter_files += (f"{toppar_dir}/{filename}",)
                    else:
                        if line.strip() != "":
                            filename = line.strip("\n").split("/")[-1]
                            parameter_files += (f"{toppar_dir}/{filename}",)

        else:
            parameter_files += (f"{toppar_dir}/top_all36_prot.rtf",)
            parameter_files += (f"{toppar_dir}/par_all36m_prot.prm",)
            parameter_files += (f"{toppar_dir}/par_all36_na.prm",)
            parameter_files += (f"{toppar_dir}/top_all36_na.rtf",)
            parameter_files += (f"{toppar_dir}/top_all36_cgenff.rtf",)
            parameter_files += (f"{toppar_dir}/par_all36_cgenff.prm",)
            parameter_files += (f"{toppar_dir}/par_all36_lipid.prm",)
            parameter_files += (f"{toppar_dir}/top_all36_lipid.rtf",)
            parameter_files += (f"{toppar_dir}/toppar_water_ions.str",)
            parameter_files += (
                f"{toppar_dir}/toppar_all36_prot_na_combined.str",
            )  # if modified aminoacids are needed
            # if os.path.isfile(f"{toppar_dir}/toppar_all36_moreions.str"):
            #     parameter_files += (f"{toppar_dir}/toppar_all36_moreions.str",)

        # set up parameter object
        parameter = pm.charmm.CharmmParameterSet(*parameter_files)
        return parameter

    def _initialize_system(
        self, configuration: dict, env: str
    ) -> pm.charmm.CharmmPsfFile:
        """
        Generates the psf file and sets the coordinates from the CHARMM-GUI files.
        Parameters
        ----------
        configuration: dict
            the configuration dictionary obtained with utils.load_config_yaml
        env: str
            waterbox,complex or vacuum
        Returns
        ----------
        psf : pm.charmm.CharmmPsfFile
        """

        if env == "vacuum":
            # take the structures from the waterbox system and extract only the ligand
            taken_from = "waterbox"
            psf_file_name = configuration["system"][self.structure][taken_from][
                "psf_file_name"
            ]
            crd_file_name = configuration["system"][self.structure][taken_from][
                "crd_file_name"
            ]

            psf_file_path = (
                f"{self.charmm_gui_base}/{taken_from}/openmm/{psf_file_name}.psf"
            )
            crd_file_path = (
                f"{self.charmm_gui_base}/{taken_from}/openmm/{crd_file_name}.crd"
            )

            # load psf
            psf = pm.charmm.CharmmPsfFile(psf_file_path)
            coord = pm.charmm.CharmmCrdFile(crd_file_path)
            psf.coordinates = coord.coordinates

            lp = False
            for atom in psf.atoms:
                if hasattr(atom, "frame_type"):
                    lp = True

            # this is used for creating the vacuum structure,
            # unfortunatly parmed forgets afterwards about the frame_type
            # which is necessary for the check_for_lp function
            if lp:
                g = psf.groups
                frame_idx = []
                frame_frame = []
                for atom in psf.atoms:
                    if hasattr(atom, "frame_type"):
                        frame_idx.append(atom.idx)
                        frame_frame.append(atom.frame_type)
                psf = psf[f":{self.tlc}"]  # the important part
                psf.groups = g
                for atom in psf.atoms:
                    if atom.idx in frame_idx:
                        atom.frame_type = frame_frame[frame_idx.index(atom.idx)]
            else:
                psf = psf[f":{self.tlc}"]

        else:
            psf_file_name = configuration["system"][self.structure][env][
                "psf_file_name"
            ]
            crd_file_name = configuration["system"][self.structure][env][
                "crd_file_name"
            ]

            psf_file_path = f"{self.charmm_gui_base}/{env}/openmm/{psf_file_name}.psf"
            crd_file_path = f"{self.charmm_gui_base}/{env}/openmm/{crd_file_name}.crd"
            psf = pm.charmm.CharmmPsfFile(psf_file_path)
            coord = pm.charmm.CharmmCrdFile(crd_file_path)
            psf.coordinates = coord.coordinates

        return psf

    def _determine_offset_and_set_possible_dummy_properties(
        self, psf: pm.charmm.CharmmPsfFile
    ) -> int:
        """
        Determines the offset and sets possible properties on the psf.
        Parameters
        ----------
        psf : pm.charmm.CharmmPsfFile
        env: str
            waterbox,complex or vacuum
        Returns
        ----------
        """
        assert type(psf) == pm.charmm.CharmmPsfFile
        if len(psf.view[f":{self.tlc}"].atoms) < 1:
            raise RuntimeError(f"No ligand selected for tlc: {self.tlc}")

        psf.number_of_dummys = 0
        psf.mutations_to_default = 0

        idx_list = []
        for atom in psf.view[f":{self.tlc}"].atoms:
            idx_list.append(int(atom.idx))

            # charge, epsilon and rmin are directly modiefied
            atom.initial_charge = atom.charge
            atom.initial_epsilon = atom.epsilon
            atom.initial_rmin = atom.rmin

        return min(idx_list)

    def _create_sdf_file(self) -> str:
        """
        Creates a sdf file, if none is available
        useful for point mutations
        """

        file_path = f"{self.charmm_gui_base}/waterbox/openmm/"

        pdb = pm.read_PDB(f"{file_path}/step3_input.pdb")

        deletedatoms = 0
        atomid = 0
        length = len(pdb.atoms)
        while deletedatoms + atomid < length:
            if pdb.atoms[atomid].residue.segid != "RNAA":
                pdb.atoms.remove(pdb.atoms[atomid])
                deletedatoms += 1
            else:
                atomid += 1

        pdb.write_pdb(f"{file_path}/step3_input_tmp.pdb")

        from openbabel import openbabel

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "sdf")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, f"{file_path}/step3_input_tmp.pdb")
        obConversion.WriteFile(mol, f"{file_path}/step3_input_reduced.sdf")

        return f"{file_path}/step3_input_reduced.sdf"

    def _return_small_molecule(self) -> Chem.rdchem.Mol:

        charmm_gui_env = self.charmm_gui_base + "waterbox"
        possible_files = []
        for ending in ["sdf", "mol", "mol2"]:
            possible_files.extend(glob.glob(f"{charmm_gui_env}/*/*{ending}"))

        # looking for small molecule
        # start with sdf
        mol2_detected = False
        mol_detected = False
        found_small_molecule = False
        for f in possible_files:
            if f.endswith(".sdf"):
                suppl = Chem.SDMolSupplier(f, removeHs=False)
                mol = next(suppl)
                logger.info(f"SDF file loaded: {f}")
                found_small_molecule = True
                return mol
            elif f.endswith(".mol2"):
                mol2_detected = True
            elif f.endswith(".mol"):
                mol_detected = True

        if mol2_detected == True or mol_detected == True:
            raise RuntimeError(
                "Please convert mol2 or mol file to sdf using obabel: {possible_files}."
            )
        if not found_small_molecule:
            raise FileNotFoundError(
                "Could not load small molecule sdf file in {charmm_gui_env}. Aborting."
            )

    def _generate_rdkit_mol(
        self, env: str, psf: pm.charmm.CharmmPsfFile
    ) -> Chem.rdchem.Mol:
        """
        Generates the rdkit mol object.
        Parameters
        ----------
        env: str
            waterbox,complex or vacuum
        psf: pm.charmm.CharmmPsfFile
        Returns
        ----------
        mol: rdkit.Chem.rdchem.Mol
        """

        assert type(psf) == pm.charmm.CharmmPsfFile

        try:
            self._create_sdf_file()
        except:
            pass

        mol = self._return_small_molecule()
        (
            atom_idx_to_atom_name,
            _,
            atom_name_to_atom_type,
            atom_idx_to_atom_partial_charge,
        ) = self.generate_atom_tables_from_psf(psf)

        for atom in mol.GetAtoms():
            atom.SetProp("atom_name", atom_idx_to_atom_name[atom.GetIdx()])
            atom.SetProp(
                "atom_type",
                atom_name_to_atom_type[atom_idx_to_atom_name[atom.GetIdx()]],
            )
            atom.SetProp("atom_index", str(atom.GetIdx()))
            atom.SetProp(
                "atom_charge", str(atom_idx_to_atom_partial_charge[atom.GetIdx()])
            )

        # check if psf and sdf have same indeces
        for a in mol.GetAtoms():
            if str(psf[a.GetIdx()].element_name) == str(a.GetSymbol()):
                pass
            else:
                raise RuntimeError("PSF to mol conversion did not work! Aborting.")

        return mol

    def generate_atom_tables_from_psf(
        self, psf: pm.charmm.CharmmPsfFile
    ) -> Tuple[dict, dict, dict, dict]:
        """
        Generate mapping dictionaries for a molecule in a psf.
        Parameters
        ----------
        psf: pm.charmm.CharmmPsfFile
        Returns
        ----------
        dict's
        """

        atom_idx_to_atom_name = dict()
        atom_name_to_atom_idx = dict()
        atom_name_to_atom_type = dict()
        atom_idx_to_atom_partial_charge = dict()

        for atom in psf.view[f":{self.tlc}"].atoms:
            atom_name = atom.name
            atom_index = atom.idx
            atom_type = atom.type
            atom_charge = atom.charge

            atom_idx_to_atom_name[atom_index] = atom_name
            atom_name_to_atom_idx[atom_name] = atom_index
            atom_name_to_atom_type[atom_name] = atom_type
            atom_idx_to_atom_partial_charge[atom_index] = atom_charge

        return (
            atom_idx_to_atom_name,
            atom_name_to_atom_idx,
            atom_name_to_atom_type,
            atom_idx_to_atom_partial_charge,
        )
