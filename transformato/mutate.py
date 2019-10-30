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
        self.s2_tlc = s2.tlc


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
            logger.info('Mol in SMILES format: {}.'.format(Chem.MolToSmiles(m,True)))

        # make copy of mols and write atom type in isotope propertie
        changed_mols = [Chem.Mol(x) for x in [m1, m2]]

        # find substructure match (ignore bond order but enforce element matching)
        mcs = rdFMCS.FindMCS(changed_mols, bondCompare=rdFMCS.BondCompare.CompareAny, timeout=120, atomCompare=rdFMCS.AtomCompare.CompareElements)
        logger.info('Substructure match: {}'.format(mcs.smartsString))

        # convert from SMARTS
        mcsp = Chem.MolFromSmarts(mcs.smartsString, False)

        s1 = (m1.GetSubstructMatch(mcsp))
        logger.info('Substructere match idx: {}'.format(s1))
        self._display_mol(m1)
        s2 = (m2.GetSubstructMatch(mcsp))
        logger.info('Substructere match idx: {}'.format(s2))
        self._display_mol(m2)

        self._substructure_match[mol1_name] = s1
        self._substructure_match[mol2_name] = s2

    def _display_mol(self, mol:Chem.Mol):
        """
        Gets mol as input and displays its 2D Structure using IPythonConsole.
        Parameters
        ----------
        mol: Chem.Mol
            a rdkit mol object
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


    def generate_mutations_to_common_core_for_mol1(self, nr_of_steps_for_el:int, nr_of_steps_for_bonded_parameters:int)->list:
        """
        Generates the mutation route to the common fore for mol1.
        Parameters
        ----------
        nr_of_steps_for_el : int
            nr of steps used for linearly scaling the charges to zero
        nr_of_steps_for_bonded_parameters : int
        Returns
        ----------
        mutations: list
            list of mutations

        """
        m = self._mutate_to_common_core('m1', self.get_common_core_idx_mol1(), nr_of_steps_for_el)
        t = self._transform_common_core(nr_of_steps_for_bonded_parameters)

        return m + t


    def generate_mutations_to_common_core_for_mol2(self, nr_of_steps_for_el:int)->list:
        """
        Generates the mutation route to the common fore for mol2.
        Returns
        ----------
        mutations: list
            list of mutations        
        """

        m = self._mutate_to_common_core('m2', self.get_common_core_idx_mol2(), nr_of_steps_for_el)
        return m

    def _transform_common_core(self, nr_of_steps_for_bonded_parameters:int)->list:
        """
        Common Core 1 is transformed to Common core 2. Bonded parameters and charges are adjusted. 
        """
        
        transformations = []
        # test if bonded mutations are necessary
        bonded_mutations_necessary = False
        for cc1, cc2 in zip(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2()):
            # did atom type change? if not don't add BondedMutations
            atom1 = self.psfs['m1'][f":{self.s1_tlc}"][cc1]
            atom2 = self.psfs['m2'][f":{self.s2_tlc}"][cc2] 
            # NOTE: This does not work!
            if atom1.type != atom2.type:
                logger.info('###########')
                logger.info('Atom that needs to be transformed: {}.'.format(atom1))
                logger.info('Atom type of atom in cc1: {}.'.format(atom1.type))
                logger.info('Template atom: {}.'.format(atom2))
                logger.info('Atom type of atom in cc2: {}.'.format(atom2.type))

                bonded_mutations_necessary = True
        
        # if necessary transform bonded parameters
        if bonded_mutations_necessary:
            logger.info('Bonded parameters need to be transformed for the cc1 topology.')
            # ugh
            t2 = BondedMutation(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2(), self.psfs['m2'], nr_of_steps_for_bonded_parameters, self.s1_tlc, self.s2_tlc)
            transformations.append(t2)
        
        # transform charges from cc1 to cc2
        logger.info('Transforming charge distribution from cc1 to cc2 in hardcoded 2 steps.')
        t1 = TransformChargesToTargetCharge(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2(), self.psfs['m2'], 2 )
        transformations.append(t1)

        return transformations


    def _mutate_to_common_core(self, name:str, cc_idx:list, nr_of_steps_for_el:int)->list:
        """
        Helper function - do not call directly.
        Generates the mutation route to the common fore for mol2.
        """

        mol = self.mols[name]
        # first LJ and electrostatic is scaled
        hydrogens = []
        mutations = []
        atoms_to_be_mutated = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in cc_idx:
                if atom.GetSymbol() == 'H':
                    hydrogens.append(idx)
                atoms_to_be_mutated.append(idx)
                logger.info('Will be decoupled: Idx:{} Element:{}'.format(idx, atom.GetSymbol()))
        
        # scale all EL of all atoms to zero
        mutations.append(ELtoZeroMutation(atoms_to_be_mutated, nr_of_steps_for_el))
        
        # scale LJ
        # start with mutation of LJ of hydrogens
        mutations.append(LJtoZeroMutation(hydrogens))
        for idx in atoms_to_be_mutated:
            if idx not in hydrogens:
                mutations.append(LJtoZeroMutation([idx]))

 
        return mutations


class BondedMutation(object):

    def __init__(self, cc1_idx:list, cc2_idx:list, cc2_psf:pm.charmm.CharmmPsfFile, nr_of_steps:int, tlc_cc1:str, tlc_cc2:str):
        """
        Scale the bonded parameters inside the common core.
        Parameters
        ----------
        cc1_idx : list
            indices of cc1
        cc2_idx : list
            indices of cc2 (in the same order as cc1)
        cc2_psf : pm.charmm.CharmmPsfFile
            the template psf that is used to generate the new bonded parmaeters
        nr_of_steps : int
        tlc_cc1 : str
            three letter code of ligand in cc1
        tlc_cc2 : str
            three letter code of ligand in cc2
        """
        self.cc1_idx = cc1_idx
        self.cc2_idx = cc2_idx
        self.cc2_psf = cc2_psf #a copy of only the ligand!
        self.nr_of_steps = nr_of_steps
        self.tlc_cc1 = tlc_cc1
        self.tlc_cc2 = tlc_cc2
        self.already_done_once = False
        self.new_atoms = []
        self.new_bonds = []
        self.new_angles = []

    def _initialize(self, psf:pm.charmm.CharmmPsfFile, offset:int):
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

        ##########################################
        # atoms matching
        for cc1, cc2 in zip(self.cc1_idx, self.cc2_idx):
            # did atom type change? if not continue
            if psf[cc1+offset].type == self.cc2_psf[cc2].type:
                continue
            atom_map_cc1_to_cc2[psf[cc1+offset].name] = self.cc2_psf[cc2].name
            psf.cc_atoms.append(psf[cc1+offset])
            self.new_atoms.append(self.cc2_psf[cc2])
            # setting real parameter 
            psf[cc1+offset].real_epsilon = psf[cc1+offset].epsilon
            psf[cc1+offset].real_sigma = psf[cc1+offset].sigma

        ##########################################
        # bonds
        for cc1_bond in psf.bonds:
            cc1_a1 = cc1_bond.atom1.name
            cc1_a2 = cc1_bond.atom2.name
            # only bonds in ligand
            if not all(elem.residue.name == self.tlc_cc1.upper() for elem in [cc1_bond.atom1, cc1_bond.atom2]):
                continue
            # all atoms of the bond must be in cc
            # everything outside the cc are bonded terms between dummies or 
            # between real atoms and dummies and we can ignore them for now
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2]):
                    continue
            # set real parameter
            cc1_bond.real_k = cc1_bond.type.k
            cc1_bond.real_req = cc1_bond.type.req

            psf.cc_bonds.append(cc1_bond)
            for cc2_bond in self.cc2_psf.bonds:
                cc2_a1 = cc2_bond.atom1.name
                cc2_a2 = cc2_bond.atom2.name
                # only bonds in ligand
                if not all(elem.residue.name == self.tlc_cc2.upper() for elem in [cc2_bond.atom1, cc2_bond.atom2]):
                    continue
                # all atoms of the bond must be in cc
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2]):
                    continue
                
                # match the two bonds
                if sorted([atom_map_cc1_to_cc2[cc1_a1], atom_map_cc1_to_cc2[cc1_a2]]) == sorted([cc2_a1, cc2_a2]):
                    self.new_bonds.append(cc2_bond)

        if(len(psf.cc_bonds) != len(self.new_bonds)):
            print('Old bonds: {}'.format(len(psf.cc_bonds)))
            print('New bonds: {}'.format(len(self.new_bonds)))
            raise RuntimeError('Nr of bonds is different in the common cores.')
        else:
            logger.info('Old bonds: {}'.format(len(psf.cc_bonds)))
            logger.info('New bonds: {}'.format(len(self.new_bonds)))

        ##########################################
        for cc1_angle in psf.angles:
            cc1_a1 = cc1_angle.atom1.name
            cc1_a2 = cc1_angle.atom2.name
            cc1_a3 = cc1_angle.atom3.name
            # only angles in ligand
            if not all(elem.residue.name == self.tlc_cc1.upper() for elem in [cc1_angle.atom1, cc1_angle.atom2, cc1_angle.atom3]):
                continue
            # only angles in cc
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2, cc1_a3]):
                    continue
            # set real parameter
            cc1_angle.real_k = cc1_angle.type.k
            cc1_angle.real_theteq = cc1_angle.type.theteq

            psf.cc_angles.append(cc1_angle)

            for cc2_angle in self.cc2_psf.angles:
                cc2_a1 = cc2_angle.atom1.name
                cc2_a2 = cc2_angle.atom2.name
                cc2_a3 = cc2_angle.atom3.name
                # only angles in ligand
                if not all(elem.residue.name == self.tlc_cc2.upper() for elem in [cc2_angle.atom1, cc2_angle.atom2, cc2_angle.atom3]):
                    continue
                # only angles in cc
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2, cc2_a3]):
                    continue

                if sorted([atom_map_cc1_to_cc2[cc1_a1], atom_map_cc1_to_cc2[cc1_a2], atom_map_cc1_to_cc2[cc1_a3]]) == sorted([cc2_a1, cc2_a2, cc2_a3]):
                    self.new_angles.append(cc2_angle)
    
        if(len(psf.cc_angles) != len(self.new_angles)):
            print('Old angles: {}'.format(len(psf.cc_angles)))
            print('New angles: {}'.format(len(self.new_angles)))
            raise RuntimeError('Nr of angles is different in the common cores.')
        else:
            logger.info('Old angles: {}'.format(len(psf.cc_angles)))
            logger.info('New angles: {}'.format(len(self.new_angles)))

        
        ##########################################
        # torsions are treated differently
        # here we add torsions from the 
        # target topology to the starting topology and 
        # turn the starting dihedrals off and the target 
        # dihedrals on

        # get all torsions present in initial topology
        for cc1_torsion in psf.dihedrals:
            cc1_a1 = cc1_torsion.atom1.name
            cc1_a2 = cc1_torsion.atom2.name
            cc1_a3 = cc1_torsion.atom3.name
            cc1_a4 = cc1_torsion.atom4.name
            # only torsions in ligand
            if not all(elem.residue.name == self.tlc_cc1.upper() for elem in [cc1_torsion.atom1, cc1_torsion.atom2, cc1_torsion.atom3, cc1_torsion.atom4]):
                continue
            # all atoms must be in the cc
            if not all(elem in atom_map_cc1_to_cc2.keys() for elem in [cc1_a1, cc1_a2, cc1_a3, cc1_a4]):
                continue
            # set identifier
            cc1_torsion.gets_modified = True
            cc1_torsion.type_cc2 = []
            for torsion_t in cc1_torsion.type:
                torsion_t.real_phi_k = torsion_t.phi_k

            # get corresponding torsion types in the new topology
            found = False
            for cc2_torsion in self.cc2_psf.dihedrals:
                cc2_a1 = cc2_torsion.atom1.name
                cc2_a2 = cc2_torsion.atom2.name
                cc2_a3 = cc2_torsion.atom3.name
                cc2_a4 = cc2_torsion.atom4.name
                # only torsions in ligand
                if not all(elem.residue.name == self.tlc_cc2.upper() for elem in [cc2_torsion.atom1, cc2_torsion.atom2, cc2_torsion.atom3, cc2_torsion.atom4]):
                    continue
                # only torsion in cc
                if not all(elem in atom_map_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
                    continue

                if sorted([atom_map_cc1_to_cc2[cc1_a1], atom_map_cc1_to_cc2[cc1_a2], atom_map_cc1_to_cc2[cc1_a3], atom_map_cc1_to_cc2[cc1_a4]]) == sorted([cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
                    found = True
                    for torsion_t in cc2_torsion.type:
                        torsion_t.real_phi_k = torsion_t.phi_k
                        cc1_torsion.type_cc2.append(torsion_t)
            if found == False:
                raise RuntimeError('Be careful - there is an unmatched torsion present.')

        ##########################################
        # impropers

        for cc1_torsion in psf.impropers:
            cc1_a1 = cc1_torsion.atom1.name
            cc1_a2 = cc1_torsion.atom2.name
            cc1_a3 = cc1_torsion.atom3.name
            cc1_a4 = cc1_torsion.atom4.name
           
            # only improper in ligand
            if not all(elem.residue.name == self.tlc_cc1.upper() for elem in [cc1_torsion.atom1, cc1_torsion.atom2, cc1_torsion.atom3, cc1_torsion.atom4]):
                continue
            # impropers don't need to have only members that are in the commen core
            # but the central atom needs to be in the commen core
            # NOTE: I think the central atom is the first! 
            if not cc1_a1 in atom_map_cc1_to_cc2.keys():
                continue
            # set real parameter
            cc1_torsion.real_psi_k = cc1_torsion.type.psi_k
            cc1_torsion.improper_present_in = 'cc1'


        for cc2_torsion in self.cc2_psf.impropers:
            cc2_a1 = cc2_torsion.atom1.name
            cc2_a2 = cc2_torsion.atom2.name
            cc2_a3 = cc2_torsion.atom3.name
            cc2_a4 = cc2_torsion.atom4.name
            # only improper in ligand
            if not all(elem.residue.name == self.tlc_cc2.upper() for elem in [cc2_torsion.atom1, cc2_torsion.atom2, cc2_torsion.atom3, cc2_torsion.atom4]):
                continue
            # impropers don't need to have only members that are in the commen core
            if not cc1_a1 in atom_map_cc1_to_cc2.keys():
                continue
            print('Found improper in cc2')
            print(cc2_torsion)
            
        #     new_improper.append(cc2_torsion)
        #     new_improper_type.append(cc2_torsion.type)

        # # append torstion and torsion type to psf
        # for tor, tor_type in zip(new_improper, new_improper_type):
        #     new_tor = deepcopy(tor)
        #     # save the initial psi_k
        #     new_tor.real_psi_k = tor.type.psi_k
        #     # set the new torsions to zero
        #     new_tor.type.psi_k = 0.0

        #     psf.impropers.append(new_tor)
        #     psf.impropers[-1].improper_present_in = 'cc2'
        #     psf.improper_types.append(deepcopy(tor_type))

        
    def mutate(self, psf:pm.charmm.CharmmPsfFile, offset:int, current_step:int):
        """
        Mutates the bonded parameters of cc1 to cc2.
        Parameters
        ----------
        psf : pm.charmm.CharmmPsfFile
            psf that gets mutated
        offset : int
            first index of ligand in the psf
        current_step : int
            the current step in the mutation protocoll
        """

        assert(type(psf) == pm.charmm.CharmmPsfFile)

        if psf.already_initialized == False:         
            self._initialize(psf, offset)
            psf.already_initialized = True

        scale = current_step/(self.nr_of_steps -1)
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
            old_angle.type.k      = (1.0 - scale) * old_angle.real_k      + scale * new_angle.type.k
            old_angle.type.theteq = (1.0 - scale) * old_angle.real_theteq + scale * new_angle.type.theteq

        # scale torsions
        for torsion in psf.dihedrals:
            # test if one of the torsions that needs to change
            if not hasattr(torsion, 'gets_modified'):
                continue
            
            # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2 
            if scale <= 0.5:
                for torsion_t in torsion.type:
                    torsion_t.phi_k = torsion_t.real_phi_k * max(((1.0 - scale * 2)), 0.0)
            else:
            # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2 
                torsion.type = deepcopy(torsion.type_cc2)
                for torsion_t in torsion.type:
                    torsion_t.phi_k = torsion_t.real_phi_k * max((scale -0.5) * 2, 0.0)

        # # scale improper
        # for torsion in psf.impropers:
        #     # test if one of the torsions that needs to change
        #     if not hasattr(torsion, 'improper_present_in'):
        #         continue
            
        #     # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2 
        #     if torsion.improper_present_in == 'cc1':
        #         # scale the initial torsions down
        #         torsion.type.psi_k = torsion.real_psi_k * max(((1.0 - scale * 2)), 0.0)
        #     # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2 
        #     elif torsion.improper_present_in == 'cc2':
        #         # scale the new torsinos up
        #         torsion.type.psi_k = torsion.real_psi_k * max((scale -0.5) * 2, 0.0)



    def _modify_type(self, atom:pm.Atom, psf:pm.charmm.CharmmPsfFile):

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
        """
        Scale the electrostatics of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_list : list
            atoms for which charges are scaled to zero
        mr_of_steps : int
            determines the scaling factor : multiplicator = 1 - (current_step / (nr_of_steps))
        """
        super().__init__(atom_idx, nr_of_steps)

    def mutate(self, psf:pm.charmm.CharmmPsfFile, offset:int, current_step:int):
        logger.info('Charges to zero mutation')
        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 1 - (current_step / (self.nr_of_steps -1))
            charge = round(atom.real_charge * multiplicator , 5)
            self._scale_charge(atom, charge)


class TransformChargesToTargetCharge(ELMutation):

    def __init__(self, atom_idx:list, target_idx:list, target_psf:pm.charmm.CharmmPsfFile, nr_of_steps:int):
        """
        Scale the electrostatics of atoms specified in the atom_idx list to the charges specified in the target_idx list.
        Parameters
        ----------
        atom_list : list
            atom indices that specify atoms for which charges are scaled
        target_idx : list
            atom indices that specify atoms that have the target charge
        target_psf : pm.charmm.CharmmPsfFile
            target_psf has the target charges 
        mr_of_steps : int
            determines the scaling factor : (1.0 - (current_step/nr_of_steps) * atom.charge + (current_step/nr_of_steps) * target_atom.charge
        """
        super().__init__(atom_idx, nr_of_steps)
        self.target_psf = target_psf
        self.target_idx = target_idx

    def mutate(self, psf, offset:int, current_step:int):
        logger.info('Mutate charge from cc1 to cc2')
        scale = current_step/(self.nr_of_steps -1)

        for idx, target_idx in zip(self.atom_idx, self.target_idx):
            atom = psf[idx + offset]
            target_atom = self.target_psf[target_idx]
            charge = (1.0 - scale) * atom.real_charge + scale * target_atom.charge
            self._scale_charge(atom, charge)



class LJtoZeroMutation(LJMutation):

    def __init__(self, atom_idx:list):
        """
        Set the VdW terms of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_list : list
            atoms for which charges are scaled to zero
        """
        super().__init__(atom_idx, 1)
    
    def mutate(self, psf, offset:int, current_step:int):
        logger.info('VdW to zero mutation')

        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 0
            self._modify_type(atom, psf)
            self._scale_epsilon(atom, multiplicator)
            self._scale_sigma(atom, multiplicator)

