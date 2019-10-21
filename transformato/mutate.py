import parmed as pm
import logging
from simtk import unit
from rdkit import Chem
import rdkit
import os, io
from copy import deepcopy   
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from IPython.core.display import display
from transformato import state

logger = logging.getLogger(__name__)


class ProposeMutationRoute(object):
    
    def __init__(self, s1:state, s2:state):
        
        """
        A class that proposes the mutation route between two molecules with a 
        common core (same atom types) based on two mols and generates the mutation 
        objects to perform the mutation on the psf objects.
        Parameters
        ----------
        mol1: Chem.Mol
        mol2: Chem.Mol
        """

        mol1_name:str = 'm1'
        mol2_name:str = 'm2'

        self.mols:dict = {mol1_name : s1.mol, mol2_name : s2.mol}
        self.psfs:dict = {mol1_name : s1.waterbox_psf[f":{s1.tlc}"], mol2_name : s2.waterbox_psf[f":{s2.tlc}"]}
        self._substructure_match:dict = { mol1_name : [], mol2_name : []}
        self._calculate_common_core(mol1_name, mol2_name)
        self.removed_indeces:dict = { mol1_name : [], mol2_name : []}
        self.added_indeces:dict = { mol1_name : [], mol2_name : []}
        self.s1_tlc = s1.tlc


    def add_common_core_atom_to_mol1(self, idx:int):
        self._add_common_core_atom('m1', idx)

    def add_common_core_atom_to_mol2(self, idx:int):
        self._add_common_core_atom('m2', idx)

    def _add_common_core_atom(self, name, idx):
        self.added_indeces[name].append(idx)
        
    def get_common_core_idx_mol1(self)->list:
        """
        Returns the common core of mol1.
        """
        return self._get_common_core('m1')
    
    def get_common_core_idx_mol2(self)->list:
        """
        Returns the common core of mol2.
        """
        return self._get_common_core('m2')

    def _get_common_core(self, name:str)->list:
        """
        Helper Function - should not be called directly.
        Returns the common core.
        """
        keep_idx = []
        for idx in set(list(self._substructure_match[name]) + self.added_indeces[name]):
            if idx not in self.removed_indeces[name]:
                keep_idx.append(idx)
        return keep_idx
                   
    def _calculate_common_core(self, mol1_name:str, mol2_name:str):
        """
        A class that proposes the mutation route between two molecules with a 
        common core (same atom types) based on two mols and generates the mutation 
        objects to perform the mutation on the psf objects.
        Parameters
        ----------
        mol1_name: str
        mol2_name: str
        """
        m1,m2 = [deepcopy(self.mols[mol1_name]), deepcopy(self.mols[mol2_name])]

        for m in [m1,m2]:
            print('Mol in SMILES format: {}.'.format(Chem.MolToSmiles(m,True)))

        # make copy of mols and write atom type in isotope propertie
        changed_mols = [Chem.Mol(x) for x in [m1, m2]]

        # find substructure match (ignore bond order but enforce element matching)
        mcs = rdFMCS.FindMCS(changed_mols, bondCompare=rdFMCS.BondCompare.CompareAny, timeout=120, atomCompare=rdFMCS.AtomCompare.CompareElements)
        print('Substructure match: {}'.format(mcs.smartsString))

        # convert from SMARTS
        mcsp = Chem.MolFromSmarts(mcs.smartsString, False)

        s1 = (m1.GetSubstructMatch(mcsp))
        print('Substructere match idx: {}'.format(s1))
        self._display_mol(m1)
        s2 = (m2.GetSubstructMatch(mcsp))
        print('Substructere match idx: {}'.format(s2))
        self._display_mol(m2)

        self._substructure_match[mol1_name] = s1
        self._substructure_match[mol2_name] = s2

    def _display_mol(self, mol:Chem.Mol):
        """
        Gets mol as input and displays its 2D Structure using IPythonConsole.
        """

        def mol_with_atom_index(mol):
            atoms = mol.GetNumAtoms()
            for idx in range(atoms):
                mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
            return mol
        
        mol = mol_with_atom_index(mol)
        AllChem.Compute2DCoords(mol)
        display(mol)


    def show_common_core_on_mol1(self):
        """
        Shows common core on mol1        
        """
        return self._show_common_core(self.mols['m1'], self.get_common_core_idx_mol1())

    def show_common_core_on_mol2(self):
        """
        Shows common core on mol2        
        """
        return self._show_common_core(self.mols['m2'], self.get_common_core_idx_mol2())

    def _show_common_core(self, mol, highlight):
        """
        Helper function - do not call directly.
        Show common core.
        """
        #https://rdkit.blogspot.com/2015/02/new-drawing-code.html
        
        mol = deepcopy(mol)
        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(600,600)
        drawer.SetFontSize(0.25)

        opts = drawer.drawOptions()

        for i in mol.GetAtoms():
            opts.atomLabels[i.GetIdx()] = str(i.GetProp('atom_index')) + ':' + i.GetProp('atom_type')

        drawer.DrawMolecule(mol, highlightAtoms=highlight)
        Draw.DrawingOptions.includeAtomNumbers=False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:','')
        return(svg)


    def generate_mutations_to_common_core_for_mol1(self)->list:
        """
        Generates the mutation route to the common fore for mol1.
        """
        m = self._mutate_to_common_core('m1', self.get_common_core_idx_mol1())
        t = [self._transform_common_core()]

        return m + t


    def generate_mutations_to_common_core_for_mol2(self)->list:
        """
        Generates the mutation route to the common fore for mol2.
        """

        m = self._mutate_to_common_core('m2', self.get_common_core_idx_mol2())
        return m

    def _transform_common_core(self):
        
        m = BondedMutation(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2(), self.psfs['m2'], 10, self.s1_tlc)
        return m


    def _mutate_to_common_core(self, name:str, cc_idx:list)->list:
        """
        Helper function - do not call directly.
        Generates the mutation route to the common fore for mol2.
        """

        mol = self.mols[name]
        # first LJ and electrostatic is scaled
        mutations = []
        atoms_to_be_mutated = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in cc_idx:
                atoms_to_be_mutated.append(idx)
                print('Will be decoupled: Idx:{} Elemeng:{}'.format(idx, atom.GetSymbol()))
        
        # scale all EL of all atoms to zero
        mutations.append(ELtoZeroMutation(atoms_to_be_mutated, 10))
        # start with mutation of VdW of hydrogens
        mutations.append(LJtoZeroMutation(atoms_to_be_mutated, 10))
 
        return mutations


class BondedMutation(object):

    def __init__(self, cc1_idx:list, cc2_idx:list, cc2_psf, nr_of_steps:int, tlc:str):
        self.cc1_idx = cc1_idx
        self.cc2_idx = cc2_idx
        self.cc2_psf = cc2_psf #a copy of only the ligand!
        self.nr_of_steps = nr_of_steps
        self.tlc = tlc
        self.already_done_once = False
        self.new_atoms = []
        self.new_bonds = []
        self.new_angles = []
        self.new_torsions = []
        self.new_improper = []


    def _initialize(self, psf, offset):
        """
        Initialize the parameter lists, also set real_{parameter} on old_atom, old_bond, ...
        """
        
        atom_map_cc1_to_cc2 = {}
        # overwritting since this is done twice #NOTE: not best coding
        self.new_atoms = []
        self.new_bonds = []
        self.new_angles = []
        self.new_torsions = []
        self.new_improper = []

        for cc1, cc2 in zip(self.cc1_idx, self.cc2_idx):
            atom_map_cc1_to_cc2[psf[cc1+offset].name] = self.cc2_psf[cc2].name
            psf.cc_atoms.append(psf[cc1+offset])
            self.new_atoms.append(self.cc2_psf[cc2])
            # setting real parameter 
            psf[cc1+offset].real_epsilon = psf[cc1+offset].epsilon
            psf[cc1+offset].real_sigma = psf[cc1+offset].sigma
        ##########################################
        for cc1_bond in psf.bonds:
            # only bonds in ligand
            if not all(elem.residue.name == self.tlc.upper() for elem in [cc1_bond.atom1, cc1_bond.atom2]):
                continue
            # only bonds in cc
            cc1_a1 = cc1_bond.atom1.name
            cc1_a2 = cc1_bond.atom2.name
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2]):
                    continue

            # set real parameter
            cc1_bond.real_k = cc1_bond.type.k
            cc1_bond.real_req = cc1_bond.type.req

            psf.cc_bonds.append(cc1_bond)
            for cc2_bond in self.cc2_psf.bonds:
                cc2_a1 = cc2_bond.atom1.name
                cc2_a2 = cc2_bond.atom2.name
                # only bonds in cc
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2]):
                    continue

                if sorted([atom_map_cc1_to_cc2[cc1_a1], atom_map_cc1_to_cc2[cc1_a2]]) == sorted([cc2_a1, cc2_a2]):
                    self.new_bonds.append(cc2_bond)

        if(len(psf.cc_bonds) != len(self.new_bonds)):
            print('Old bonds: {}'.format(len(psf.cc_bonds)))
            print('New bonds: {}'.format(len(self.new_bonds)))
            raise RuntimeError('Nr of bonds is different in the common cores.')
        
        ##########################################
        for cc1_angle in psf.angles:
            cc1_a1 = cc1_angle.atom1.name
            cc1_a2 = cc1_angle.atom2.name
            cc1_a3 = cc1_angle.atom3.name
            # only angles in ligand
            if not all(elem.residue.name == self.tlc.upper() for elem in [cc1_angle.atom1, cc1_angle.atom2, cc1_angle.atom3]):
                continue
            # only angles in cc
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2, cc1_a3]):
                    continue

            psf.cc_angles.append(cc1_angle)

            for cc2_angle in self.cc2_psf.angles:
                cc2_a1 = cc2_angle.atom1.name
                cc2_a2 = cc2_angle.atom2.name
                cc2_a3 = cc2_angle.atom3.name
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2, cc2_a3]):
                    continue

                if sorted([atom_map_cc1_to_cc2[cc1_a1], atom_map_cc1_to_cc2[cc1_a2], atom_map_cc1_to_cc2[cc1_a3]]) == sorted([cc2_a1, cc2_a2, cc2_a3]):
                    self.new_angles.append(cc2_angle)
    
        if(len(psf.cc_angles) != len(self.new_angles)):
            raise RuntimeError('Nr of angles is different in the common cores.')

        ##########################################
        for cc1_torsion in psf.dihedrals:
            cc1_a1 = cc1_torsion.atom1.name
            cc1_a2 = cc1_torsion.atom2.name
            cc1_a3 = cc1_torsion.atom3.name
            cc1_a4 = cc1_torsion.atom4.name
            # only torsions in ligand
            if not all(elem.residue.name == self.tlc.upper() for elem in [cc1_torsion.atom1, cc1_torsion.atom2, cc1_torsion.atom3, cc1_torsion.atom4]):
                continue
            # only torsions in cc
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2, cc1_a3, cc1_a4]):
                continue

            psf.cc_torsions.append(cc1_torsion)

            for cc2_torsion in self.cc2_psf.dihedrals:
                cc2_a1 = cc2_torsion.atom1.name
                cc2_a2 = cc2_torsion.atom2.name
                cc2_a3 = cc2_torsion.atom3.name
                cc2_a4 = cc2_torsion.atom4.name
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
                    continue
                if sorted([atom_map_cc1_to_cc2[cc1_a1], atom_map_cc1_to_cc2[cc1_a2], atom_map_cc1_to_cc2[cc1_a3], atom_map_cc1_to_cc2[cc1_a4]]) == sorted([cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
                    self.new_torsions.append(cc2_torsion)

        if(len(psf.cc_torsions) != len(self.new_torsions)):
            raise RuntimeError('Nr of torsions is different in the common cores: {} = {}.'.format(len(torsion_list_cc1), len(torsion_list_cc2)))

        ##########################################
        for cc1_torsion in psf.impropers:
            cc1_a1 = cc1_torsion.atom1.name
            cc1_a2 = cc1_torsion.atom2.name
            cc1_a3 = cc1_torsion.atom3.name
            cc1_a4 = cc1_torsion.atom4.name
                # only torsions in ligand
            if not all(elem.residue.name == self.tlc.upper() for elem in [cc1_torsion.atom1, cc1_torsion.atom2, cc1_torsion.atom3, cc1_torsion.atom4]):
                continue
            # only torsions in cc
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2, cc1_a3, cc1_a4]):
                continue

            psf.cc_improper.append(cc1_torsion)

            for cc2_torsion in self.cc2_psf.impropers:
                cc2_a1 = cc2_torsion.atom1.name
                cc2_a2 = cc2_torsion.atom2.name
                cc2_a3 = cc2_torsion.atom3.name
                cc2_a4 = cc2_torsion.atom4.name
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
                    continue
                psf.cc_improper.append(cc2_torsion)

        
    def mutate(self, psf, offset, current_step:int):

        if psf.already_done_once == False:         
            self._initialize(psf, offset)
            psf.already_done_once = True

        scale = current_step/(self.nr_of_steps - 1)
        print('Bonded parameters scaled.')
        # scale atoms
        for old_atom, new_atom in zip(psf.cc_atoms, self.new_atoms):
            self._modify_type(old_atom, psf)
            old_atom.epsilon = (1.0 - scale) * old_atom.real_epsilon + scale * new_atom.epsilon
            old_atom.sigma = (1.0 - scale) * old_atom.real_sigma + scale * new_atom.sigma

        # scale bonds
        for old_bond, new_bond in zip(psf.cc_bonds, self.new_bonds):
            old_bond.type.k   = (1.0 - scale) * old_bond.real_k   + scale * new_bond.type.k
            old_bond.type.req = (1.0 - scale) * old_bond.real_req + scale * new_bond.type.req

        # scale angles
        for old_angle, new_angle in zip(psf.cc_angles, self.new_angles):
            old_angle.type.k      = (1.0 - scale) * old_angle.type.k      + scale * new_angle.type.k
            old_angle.type.theteq = (1.0 - scale) * old_angle.type.theteq + scale * new_angle.type.theteq

        # scale torsions
        # for old_torsion, new_torsion in zip(psf.cc_torsions, self.new_torsions):
        #      if current_step < (self.nr_of_step/2):
        #          for t in old_torsion:


            
        #     old_torsion.type.k      = (1.0 - scale) * old_torsion.type.k    + scale * new_torsion.type.k
        #     old_torsion.type.theteq = (1.0 - scale) * old_angle.type.theteq + scale * new_torsion.type.theteq




    def _modify_type(self, atom, psf):

        if (hasattr(atom, 'real_type')):
            # only change parameters
            pass
        else:
            print('Setting RRR atomtype.')
            atom.real_type = atom.type
            atom.type = f"RRR{psf.number_of_dummys}"
            psf.number_of_dummys += 1





class BaseMutation(object):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        assert(type(atom_idx) == list)
        self.atom_idx = atom_idx
        self.nr_of_steps = nr_of_steps

class ELMutation(BaseMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)

    def _scale_charge(self, atom, charge):
        atom.charge = charge


class LJMutation(BaseMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)

    def _scale_epsilon(self, atom, multiplicator):
        epsilon = atom.real_epsilon * multiplicator
        atom.epsilon = epsilon

    def _scale_sigma(self, atom, multiplicator):
        sigma = atom.real_sigma *  multiplicator
        atom.sigma = sigma

    def _modify_type(self, atom, psf):

        if (hasattr(atom, 'real_type')):
            # only change parameters
            pass
        else:
            atom.real_type = atom.type
            atom.type = f"DDD{psf.number_of_dummys}"
            psf.number_of_dummys += 1



class ELtoZeroMutation(ELMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)

    def mutate(self, psf, offset:int, current_step:int):
        print('EL to Zero Mutation')
        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 1 - (current_step / (self.nr_of_steps -1))
            charge = round(atom.real_charge * multiplicator , 5)
            self._scale_charge(atom, charge)


class LJtoZeroMutation(LJMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)
    
    def mutate(self, psf, offset:int, current_step:int):
        print('VdW to Zero Mutation')

        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 0
            self._modify_type(atom, psf)
            self._scale_epsilon(atom, multiplicator)
            self._scale_sigma(atom, multiplicator)



class LJtoZeroMutationWithSoftCore(LJMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)
    
    def mutate(self, psf, offset:int, current_step:int):
        print('VdW to Zero Mutation')

        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 1 - (current_step / (self.nr_of_steps -1))
            self._modify_type(atom, psf)
            self._scale_epsilon(atom, multiplicator)
            self._scale_sigma(atom, multiplicator)



