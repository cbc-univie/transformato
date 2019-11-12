import parmed as pm
import logging
from simtk import unit
from rdkit import Chem
import rdkit
import os, io
from .utils import get_toppar_dir
from collections import namedtuple

logger = logging.getLogger(__name__)


class SystemStructure(object):

    def __init__(self, configuration:dict, structure:str):
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


        self.structure:str = structure
        self.name:str = configuration['system'][structure]['name']
        self.tlc:str = configuration['system'][structure]['tlc']
        self.charmm_gui_base:str = configuration['system'][structure]['charmm_gui_dir']

        # running a binding-free energy calculation? 
        if configuration['simulation']['free-energy-type'] == 'binding-free-energy':
            self.envs:list = ['complex', 'waterbox']
            self.parameter:pm.charmm.CharmmParameterSet = self._read_parameters(configuration, 'complex')

            # set up complex objects
            self.complex_psf:pm.charmm.CharmmPsfFile = self._initialize_system(configuration, 'complex')
            # load parameters
            self.complex_psf.load_parameters(self.parameter)
            # get offset
            self.complex_offset:int = self._determine_offset_and_set_possible_dummy_properties(self.complex_psf)
            
            # set up waterbox objects
            self.waterbox_psf:pm.charmm.CharmmPsfFile = self._initialize_system(configuration, 'waterbox')
            # load parameters
            self.waterbox_psf.load_parameters(self.parameter)
            # get offset
            self.waterbox_offset:int = self._determine_offset_and_set_possible_dummy_properties(self.waterbox_psf)

            # generate rdkit mol object of small molecule
            self.mol:Chem.Mol = self._generate_rdkit_mol('complex', self.complex_psf[f":{self.tlc}"])

        else:
            raise NotImplementedError('solvation free energy not finished yet.')
            self.envs = ['waterbox', 'vacuum']
            self.parameter = self._read_parameters(configuration, 'waterbox')


    def _read_parameters(self, configuration:dict, env:str)->pm.charmm.CharmmParameterSet:
        """
        Reads in topparameters from a toppar dir and ligand specific parameters.
        Parameters
        ----------
        configuration: dict
            the configuration dictionary obtained with utils.load_config_yaml
        env: str
            waterbox,complex or vacuum
        Returns
        ----------
        parameter : pm.charmm.CharmmParameterSet
            parameters obtained from the CHARMM-GUI output dir.      
        """
 
        charmm_gui_env = self.charmm_gui_base + env
        tlc = self.tlc       
        tlc_lower = str(tlc).lower()
        parameter_files = f"{charmm_gui_env}/{tlc_lower}/{tlc_lower}.rtf", f"{charmm_gui_env}/{tlc_lower}/{tlc_lower}.prm"         
        toppar_dir = get_toppar_dir()
        parameter_files =  parameter_files + (f"{toppar_dir}/top_all36_cgenff.rtf",) 
        parameter_files =  parameter_files + (f"{toppar_dir}/par_all36_cgenff.prm",)
        parameter_files =  parameter_files + (f"{toppar_dir}/top_all36_prot.rtf",)
        parameter_files =  parameter_files + (f"{toppar_dir}/par_all36m_prot.prm",)
        parameter_files =  parameter_files + (f"{toppar_dir}/par_all36_na.prm",)
        parameter_files =  parameter_files + (f"{toppar_dir}/top_all36_na.rtf",)
        parameter_files =  parameter_files + (f"{toppar_dir}/par_all36_lipid.prm",)
        parameter_files =  parameter_files + (f"{toppar_dir}/par_all36_na.prm",)
        parameter_files =  parameter_files + (f"{toppar_dir}/top_all36_lipid.rtf",)
        parameter_files =  parameter_files + (f"{toppar_dir}/toppar_water_ions.str",)
            
        # set up parameter objec
        parameter = pm.charmm.CharmmParameterSet(*parameter_files)
        return parameter

    def _initialize_system(self, configuration:dict, env:str)->pm.charmm.CharmmPsfFile:
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
        
        psf_file_name = configuration['system'][self.structure][env]['psf_file_name']
        crd_file_name = configuration['system'][self.structure][env]['crd_file_name']

        psf_file_path = f"{self.charmm_gui_base}/{env}/openmm/{psf_file_name}.psf"
        crd_file_path = f"{self.charmm_gui_base}/{env}/openmm/{crd_file_name}.crd"
        psf = pm.charmm.CharmmPsfFile(psf_file_path)
        coord = pm.charmm.CharmmCrdFile(crd_file_path)
        psf.coordinates = coord.coordinates

        return psf


    def _determine_offset_and_set_possible_dummy_properties(self, psf:pm.charmm.CharmmPsfFile)->int:
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
        assert(type(psf) == pm.charmm.CharmmPsfFile)
        psf.number_of_dummys = 0

        idx_list = []
        for atom in psf.view[f":{self.tlc}"].atoms:
            idx_list.append(int(atom.idx))
            
            # charge, epsilon and rmin are directly modiefied
            atom.initial_charge = atom.charge
            atom.initial_epsilon = atom.epsilon
            atom.initial_rmin = atom.rmin

        return min(idx_list)
    
    def _generate_rdkit_mol(self, env:str, psf:pm.charmm.CharmmPsfFile)->Chem.Mol:
        """
        Generates the rdkit mol object.
        Parameters
        ----------
        env: str
            waterbox,complex or vacuum
        psf: pm.charmm.CharmmPsfFile
        Returns
        ----------
        mol: rdkit.Chem.mol
        """

        assert(type(psf) == pm.charmm.CharmmPsfFile)
        charmm_gui_env = self.charmm_gui_base + env
        tlc = self.tlc       
        tlc_lower = str(tlc).lower()
        sdf_file = f"{charmm_gui_env}/{tlc_lower}/{tlc_lower}.sdf"        

        
        mol = Chem.MolFromMolFile(sdf_file, removeHs=False)
        atom_idx_to_atom_name, _, atom_name_to_atom_type, atom_idx_to_atom_partial_charge = self.generate_atom_tables_from_psf(psf)

        for atom in mol.GetAtoms():
            atom.SetProp('atom_name', atom_idx_to_atom_name[atom.GetIdx()])
            atom.SetProp('atom_type', atom_name_to_atom_type[atom_idx_to_atom_name[atom.GetIdx()]])
            atom.SetProp('atom_index', str(atom.GetIdx()))
            atom.SetProp('atom_charge', str(atom_idx_to_atom_partial_charge[atom.GetIdx()]))
            
        # check if psf and sdf have same indeces
        for a in mol.GetAtoms():
            if str(psf[a.GetIdx()].element_name) == str(a.GetSymbol()):
                pass
            else:
                raise RuntimeError('PSF to mol conversion did not work! Aborting.')

        return mol


    def generate_atom_tables_from_psf(self, psf:pm.charmm.CharmmPsfFile)->(dict,dict,dict,dict):
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

        for atom in psf.view[':' + str(self.tlc)].atoms:
            atom_name = atom.name
            atom_index = atom.idx
            atom_type = atom.type
            atom_charge = atom.charge

            atom_idx_to_atom_name[atom_index] = atom_name
            atom_name_to_atom_idx[atom_name] = atom_index
            atom_name_to_atom_type[atom_name] = atom_type
            atom_idx_to_atom_partial_charge[atom_index] = atom_charge

        return (atom_idx_to_atom_name, atom_name_to_atom_idx, atom_name_to_atom_type, atom_idx_to_atom_partial_charge)


