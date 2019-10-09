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

logger = logging.getLogger(__name__)


class ProposeMutationRoute(object):
    
    def __init__(self, mol1:Chem.Mol, mol2:Chem.Mol):
        
        """
        A class that proposes the mutation route between two molecules with a 
        commen core (same atom types) based on two mols and generates the mutation 
        objects to perform the mutation on the psf objects.
        Parameters
        ----------
        mol1: Chem.Mol
        mol2: Chem.Mol
        """

        mol1_name:str = 'm1'
        mol2_name:str = 'm2'

        self.mols:dict = {mol1_name : mol1, mol2_name : mol2}
        self._substructure_match:dict = { mol1_name : [], mol2_name : []}
        self._calculate_commen_core(mol1_name, mol2_name)
        self.removed_indeces:dict = { mol1_name : [], mol2_name : []}
        self.added_indeces:dict = { mol1_name : [], mol2_name : []}

    def get_commen_core_idx_mol1(self)->list:
        """
        Returns the commen core of mol1.
        """
        return self._get_commen_core('m1')
    
    def get_commen_core_idx_mol2(self)->list:
        """
        Returns the commen core of mol2.
        """
        return self._get_commen_core('m2')

    def _get_commen_core(self, name:str)->list:
        """
        Helper Function - should not be called directly.
        Returns the commen core of mol2.
        """
        keep_idx = []
        for idx in set(list(self._substructure_match[name]) + self.added_indeces[name]):
            if idx not in self.removed_indeces[name]:
                keep_idx.append(idx)
        return keep_idx
                   
    def _calculate_commen_core(self, mol1_name:str, mol2_name:str):
        """
        A class that proposes the mutation route between two molecules with a 
        commen core (same atom types) based on two mols and generates the mutation 
        objects to perform the mutation on the psf objects.
        Parameters
        ----------
        mol1_name: str
        mol2_name: str
        """
        m1,m2 = [deepcopy(self.mols[mol1_name]), deepcopy(self.mols[mol2_name])]
        atom_type_hash = dict()
        replacment_dict = dict()
        replacment_dict['*'] = '' 

        for m in [m1,m2]:
            print('Mol in SMILES format: {}.'.format(Chem.MolToSmiles(m,True)))

        # make copy of mols and write atom type in isotope propertie
        changed_mols = [Chem.Mol(x) for x in [m1, m2]]

        for m in changed_mols:
            for idx, atom in enumerate(m.GetAtoms()):
                if atom.GetProp('atom_type') in atom_type_hash:
                    atom.SetIsotope(atom_type_hash[atom.GetProp('atom_type')])
                else:
                    isotope_value = 100 + idx
                    atom_type_hash[atom.GetProp('atom_type')] = isotope_value

                    replacment_dict[str(isotope_value)] = atom.GetSymbol() 
                    atom.SetIsotope(atom_type_hash[atom.GetProp('atom_type')])

        # find substructure match (ignore bond order but enforce atom type matching)
        mcs = rdFMCS.FindMCS(changed_mols, bondCompare=rdFMCS.BondCompare.CompareAny, timeout=120, atomCompare=rdFMCS.AtomCompare.CompareIsotopes)
        print('Substructure match: {}'.format(mcs.smartsString))

        # do some conversion between SMARTS and SMILES to get SubstructureMatching working
        # also reconvert the SMARTS pattern using the replacment_dict conversinos
        mcsp = Chem.MolFromSmarts(mcs.smartsString, False, replacements = replacment_dict)
        g = Chem.MolFromSmiles(Chem.MolToSmiles(mcsp, allHsExplicit=True), sanitize=False)
        print('Substructure match: {}'.format(Chem.MolToSmiles(g)))

        s1 = (m1.GetSubstructMatch(g))
        self._display_mol(m1)
        s2 = (m2.GetSubstructMatch(g))
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


    def show_commen_core_on_mol1(self):
        """
        Shows commen core on mol1        
        """
        return self._show_commen_core(self.mols['m1'], self.get_commen_core_idx_mol1())

    def show_commen_core_on_mol2(self):
        """
        Shows commen core on mol2        
        """
        return self._show_commen_core(self.mols['m2'], self.get_commen_core_idx_mol2())

    def _show_commen_core(self, mol, highlight):
        """
        Helper function - do not call directly.
        Show commen core.
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


    def generate_mutations_to_commen_core_for_mol1(self)->list:
        """
        Generates the mutation route to the commen fore for mol1.
        """
        return self._mutate_to_commen_core(self.mols['m1'], self.get_commen_core_idx_mol1())


    def generate_mutations_to_commen_core_for_mol2(self)->list:
        """
        Generates the mutation route to the commen fore for mol2.
        """

        return self._mutate_to_commen_core(self.mols['m2'], self.get_commen_core_idx_mol2())

    def _mutate_to_commen_core(self, mol:Chem.Mol, cc_idx:list)->list:
        """
        Helper function - do not call directly.
        Generates the mutation route to the commen fore for mol2.
        """

        mutations = []

        atoms_to_be_mutated = []
        hydrogens = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in cc_idx:
                atoms_to_be_mutated.append(idx)
                print(atom.GetSymbol())
                if atom.GetSymbol() == 'H':
                    hydrogens.append(idx)
                print('Needs to be mutated: ', idx)
        
        # scale all EL of all atoms to zero
        mutations.append(ELtoZeroMutation(atoms_to_be_mutated, 10))
        # start with mutation of VdW of hydrogens
        for h in hydrogens:
            mutations.append(VdWtoZeroMutation([h], 1))
        # continue with heavy atoms
        for a in atoms_to_be_mutated:
            if a not in hydrogens:
                mutations.append(VdWtoZeroMutation([a], 1))
 
        return mutations


    def _mutate_cc(self, mol, cc_idx):
        pass
   



class BaseMutation(object):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        assert(type(atom_idx) == list)
        self.atom_idx = atom_idx
        self.nr_of_steps = nr_of_steps

class ELMutation(BaseMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)

    def _scale_charge(self, atom, charge):
    
        #print('Old charge: {}'.format(old_charge))
        atom.charge = charge
        #print('New charge: {}'.format(new_charge))


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


class VdWMutation(BaseMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)

    def _scale_epsilon(self, atom, epsilon):
        atom.epsilon = epsilon

    def _scale_sigma(self, atom, sigma):
        atom.sigma = sigma


    def _modify_type(self, atom, psf):

        if (hasattr(atom, 'real_type')):
            # only change parameters
            pass
        else:
            atom.real_type = atom.type
            atom.type = f"DDD{psf.number_of_dummys}"
            psf.number_of_dummys += 1

class VdWtoZeroMutation(VdWMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)
    
    def mutate(self, psf, offset:int, current_step:int):
        print('VdW to Zero Mutation')

        for i in self.atom_idx:
            atom = psf[i + offset]
            self._modify_type(atom, psf)
            epsilon = atom.real_epsilon * 0.0
            sigma = atom.real_sigma *  0.0
            self._scale_epsilon(atom, epsilon)
            self._scale_sigma(atom, sigma)





class BondedMutation(object):

    def __init__(self, atom_idx:int, nr_of_steps):

        self.atom_idx = atom_idx
        self.nr_of_steps = nr_of_steps

    # def mutate():
    #     for bond in psf[':' + tlc].bonds:
    #         bond.changed_p = ['...']




