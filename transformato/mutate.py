import io
import logging
import os
from collections import namedtuple
from copy import deepcopy

import numpy as np
import parmed as pm
import rdkit
from IPython.core.display import display
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
from simtk import unit
from transformato.system import SystemStructure
import networkx as nx 
from transformato.state import IntermediateStateFactory

logger = logging.getLogger(__name__)


class ProposeMutationRoute(object):

    def __init__(self, s1: SystemStructure, s2: SystemStructure):
        """
        A class that proposes the mutation route between two molecules with a 
        common core (same atom types) based on two mols and generates the mutation 
        objects to perform the mutation on the psf objects.
        Parameters
        ----------
        mol1: Chem.Mol
        mol2: Chem.Mol
        """

        mol1_name: str = 'm1'
        mol2_name: str = 'm2'
        self.system: dict = {'system1' : s1, 'system2': s2}
        self.mols: dict = {mol1_name: s1.mol, mol2_name: s2.mol}
        self.graphs: dict = {mol1_name: s1.graph, mol2_name: s2.graph}
        self.psfs: dict = {mol1_name: s1.waterbox_psf[f":{s1.tlc}"], mol2_name: s2.waterbox_psf[f":{s2.tlc}"]}
        self._substructure_match: dict = {mol1_name: [], mol2_name: []}
        self._calculate_common_core(mol1_name, mol2_name)
        self.removed_indeces: dict = {mol1_name: [], mol2_name: []}
        self.added_indeces: dict = {mol1_name: [], mol2_name: []}
        self.s1_tlc = s1.tlc
        self.s2_tlc = s2.tlc

        terminal_atoms_cc1, real_atoms_cc1 = self._find_terminal_atom(self.get_common_core_idx_mol1(), self.mols['m1'])
        terminal_atoms_cc2, real_atoms_cc2 = self._find_terminal_atom(self.get_common_core_idx_mol2(), self.mols['m2'])

        self.real_atom_cc1 = real_atoms_cc1[0]
        self.real_atom_cc2 = real_atoms_cc2[0]       
        self.terminal_atom_cc1 = terminal_atoms_cc1[0]
        self.terminal_atom_cc2 = terminal_atoms_cc2[0]
        self.charge_compensated_ligand1_psf, self.charge_compensated_ligand2_psf = self._prepare_cc_for_charge_transfer()


    def _redo(self):
        # these need to be recomputed when changing the common core idxs
        
        terminal_atoms_cc1, real_atoms_cc1 = self._find_terminal_atom(self.get_common_core_idx_mol1(), self.mols['m1'])
        terminal_atoms_cc2, real_atoms_cc2 = self._find_terminal_atom(self.get_common_core_idx_mol2(), self.mols['m2'])

        self.real_atom_cc1 = real_atoms_cc1[0]
        self.real_atom_cc2 = real_atoms_cc2[0]       
        self.terminal_atom_cc1 = terminal_atoms_cc1[0]
        self.terminal_atom_cc2 = terminal_atoms_cc2[0]
        self.charge_compensated_ligand1_psf, self.charge_compensated_ligand2_psf = self._prepare_cc_for_charge_transfer()

    def _prepare_cc_for_charge_transfer(self):
        # we have to run the same charge mutation that will be run on cc2 to get the 
        # charge distribution AFTER the full mutation

        # mare a copy of the psf
        m2_psf = self.psfs['m2'][:,:,:]
        m1_psf = self.psfs['m1'][:,:,:]
        charge_transformed_psfs = []
        for psf, tlc, cc_idx, real_atom in zip([m1_psf, m2_psf], 
                                               [self.s1_tlc, self.s2_tlc], 
                                               [self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2()],
                                               [self.real_atom_cc1, self.real_atom_cc2]):
            # set `initial_charge` parameter for ChargeToZeroMutation
            for atom in psf.view[f":{tlc.upper()}"].atoms:           
                # charge, epsilon and rmin are directly modiefied
                atom.initial_charge = atom.charge

            offset = min([atom.idx for atom in psf.view[f":{tlc.upper()}"].atoms])
            
            # getting copy of the atoms
            atoms_to_be_mutated = []
            for atom in psf.view[f":{tlc.upper()}"].atoms:
                idx = atom.idx - offset
                if idx not in cc_idx:
                    atoms_to_be_mutated.append(idx)
            logger.info('############################')
            logger.info('Preparing cc2 for charge transfer')
            logger.info(f"Atoms for which charge is set to zero: {atoms_to_be_mutated}")
            m = ChargeToZeroMutation(atoms_to_be_mutated, 1, cc_idx, real_atom)
            m.mutate(psf, tlc, 1)
            charge_transformed_psfs.append(psf)
        return charge_transformed_psfs[0], charge_transformed_psfs[1]


    def generate_mutation_list(self):
        
        # there are three obvious cases that we want to distinquish:
        # 1) mol1 is in mol2 (Methane -- Ethane)
        
        mutation_list = self.generate_mutations_to_common_core_for_mol1(
            nr_of_steps_for_el=5, nr_of_steps_for_cc_transformation=2)
        # write intermediate states for systems
        intermediate_state = IntermediateStateFactory(system=self.system['system1'], mutation_list=mutation_list, configuration=configuration)
        intermediate_state.generate_intermediate_states()
        
        # generate mutation route
        mutation_list = self.generate_mutations_to_common_core_for_mol2(nr_of_steps_for_el=5)
        # write intermediate states
        intermediate_state = IntermediateStateFactory(system=self.system['system2'], mutation_list=mutation_list, configuration=configuration)
        intermediate_state.generate_intermediate_states()
    
    def remove_idx_from_common_core_of_mol1(self, idx:int):
        self._remove_idx_from_common_core('m1', idx)
    

    def remove_idx_from_common_core_of_mol2(self, idx:int):
        self._remove_idx_from_common_core('m2', idx)


    def _remove_idx_from_common_core(self, name:str, idx: int):
        if idx in self.added_indeces[name] or idx in self._get_common_core(name):
            self.removed_indeces[name].append(idx)
        else:
            print(f"Idx: {idx} not in common core.")

    def add_idx_to_common_core_of_mol1(self, idx: int):
        self._add_common_core_atom('m1', idx)
        self._redo()
        print(self.get_common_core_idx_mol1())

    def add_idx_to_common_core_of_mol2(self, idx: int):
        self._add_common_core_atom('m2', idx)
        self._redo()
        print(self.get_common_core_idx_mol2())

    def _add_common_core_atom(self, name: str , idx: int):
        if idx in self.added_indeces[name] or idx in self._get_common_core(name):
            print(f"Idx: {idx} already in common core.")
            pass
        self.added_indeces[name].append(idx)

    def get_common_core_idx_mol1(self) -> list:
        """
        Returns the common core of mol1.
        """
        return self._get_common_core('m1')

    def get_common_core_idx_mol2(self) -> list:
        """
        Returns the common core of mol2.
        """
        return self._get_common_core('m2')

    def _get_common_core(self, name: str) -> list:
        """
        Helper Function - should not be called directly.
        Returns the common core.
        """
        keep_idx = []
        for idx in self._substructure_match[name] + self.added_indeces[name]:  # BEWARE: the ordering is important - don't cast set!
            if idx not in self.removed_indeces[name]:
                keep_idx.append(idx)
        return keep_idx

    def _calculate_common_core(self, mol1_name: str, mol2_name: str):
        """
        A class that proposes the mutation route between two molecules with a 
        common core (same atom types) based on two mols and generates the mutation 
        objects to perform the mutation on the psf objects.
        Parameters
        ----------
        mol1_name: str
        mol2_name: str
        """
        m1, m2 = [deepcopy(self.mols[mol1_name]), deepcopy(self.mols[mol2_name])]

        for m in [m1, m2]:
            logger.info('Mol in SMILES format: {}.'.format(Chem.MolToSmiles(m, True)))

        # make copy of mols
        changed_mols = [Chem.Mol(x) for x in [m1, m2]]

        # find substructure match (ignore bond order but enforce element matching)
        mcs = rdFMCS.FindMCS(changed_mols, bondCompare=rdFMCS.BondCompare.CompareOrder,
                             timeout=120, atomCompare=rdFMCS.AtomCompare.CompareElements)
        logger.info('Substructure match: {}'.format(mcs.smartsString))

        # convert from SMARTS
        mcsp = Chem.MolFromSmarts(mcs.smartsString, False)

        s1 = (m1.GetSubstructMatch(mcsp))
        logger.info('Substructere match idx: {}'.format(s1))
        self._display_mol(m1)
        s2 = (m2.GetSubstructMatch(mcsp))
        logger.info('Substructere match idx: {}'.format(s2))
        self._display_mol(m2)

        self._substructure_match[mol1_name] = list(s1)
        self._substructure_match[mol2_name] = list(s2)

    def _display_mol(self, mol: Chem.Mol):
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
        # https://rdkit.blogspot.com/2015/02/new-drawing-code.html

        mol = deepcopy(mol)
        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(800, 800)
        drawer.SetFontSize(0.3)

        opts = drawer.drawOptions()

        for i in mol.GetAtoms():
            opts.atomLabels[i.GetIdx()] = str(i.GetProp('atom_index')) + ':' + i.GetProp('atom_type')

        drawer.DrawMolecule(mol, highlightAtoms=highlight)
        Draw.DrawingOptions.includeAtomNumbers = False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:', '')
        return(svg)

    def generate_mutations_to_common_core_for_mol1(self, nr_of_steps_for_el: int, nr_of_steps_for_cc_transformation: int) -> list:
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
        t = self._transform_common_core(nr_of_steps_for_cc_transformation)

        return m + t

    def generate_mutations_to_common_core_for_mol2(self, nr_of_steps_for_el: int) -> list:
        """
        Generates the mutation route to the common fore for mol2.
        Returns
        ----------
        mutations: list
            list of mutations        
        """

        m = self._mutate_to_common_core('m2', self.get_common_core_idx_mol2(), nr_of_steps_for_el)
        return m

    def _transform_common_core(self, nr_of_steps_for_cc_transformation: int) -> list:
        """
        Common Core 1 is transformed to Common core 2. Bonded parameters and charges are adjusted. 
        """

        transformations = []
        logger.info('##############################')
        logger.info('##############################')
        logger.info('Transform common core')
        logger.info('##############################')
        logger.info('##############################')

        # test if bonded mutations are necessary
        bonded_terms_mutation = False
        charge_mutation = False
        for cc1, cc2 in zip(self.get_common_core_idx_mol1() + [self.terminal_atom_cc1], self.get_common_core_idx_mol2() + [self.terminal_atom_cc2]):
            # did atom type change? if not don't add BondedMutations
            atom1 = self.psfs['m1'][cc1]
            print(atom1, atom1.type)
            atom2 = self.psfs['m2'][cc2]
            print(atom2, atom2.type)
            if atom1.type != atom2.type:
                logger.info('##############################')
                logger.info('Atom type transformation')
                logger.info(f'Atom that needs to be transformed: {atom1}.')
                logger.info(f'Atom type of atom in cc1: {atom1.type}.')
                logger.info(f'Template atom: {atom2}.')
                logger.info(f'Atom type of atom in cc2: {atom2.type}.')
                bonded_terms_mutation = True
        
        for cc1, cc2 in zip(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2()):
            atom1 = self.charge_compensated_ligand1_psf[cc1]
            atom2 = self.charge_compensated_ligand2_psf[cc2]
            if atom1.charge != atom2.charge:
                logger.info('##############################')
                logger.info('Charge transformation')
                logger.info('Charge needs to be transformed on common core')
                logger.info(f'Atom that needs to be transformed: {atom1}.')
                logger.info(f'Atom charge of atom in cc1: {atom1.charge}.')
                logger.info(f'Template atom: {atom2}.')
                logger.info(f'Atom charge of atom in cc2: {atom2.charge}.')
                charge_mutation = True 
        
        
        # if necessary transform bonded parameters
        if bonded_terms_mutation or charge_mutation:
            logger.info(f'Bonded parameters mutation: {bonded_terms_mutation}.')
            logger.info(f'Charge parameters mutation: {charge_mutation}.')

            t = CommonCoreTransformation(
                self.get_common_core_idx_mol1(),
                self.get_common_core_idx_mol2(), 
                self.psfs['m1'], 
                self.psfs['m2'],
                nr_of_steps_for_cc_transformation, 
                self.s1_tlc, 
                self.s2_tlc, 
                self.terminal_atom_cc1,
                self.terminal_atom_cc2,
                self.charge_compensated_ligand2_psf[cc2],
                charge_mutation=charge_mutation,
                bonded_terms_mutation=bonded_terms_mutation)
            transformations.append(t)
        else:
            logger.info('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            logger.info('No transformations needed.')
            logger.info('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            transformations = []

        return transformations

    @staticmethod
    def _find_terminal_atom(cc_idx: list, mol:Chem.Mol):
        """
        Find atoms that connect the rest of the molecule to the common core.

        Args:
            cc_idx (list): common core index atoms
            mol ([type]): rdkit mol object
        """
        terminal_atoms = []
        last_real_atoms = []
        
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in cc_idx:
                neighbors = [x.GetIdx() for x in atom.GetNeighbors()]
                if any([n in cc_idx for n in neighbors]):
                    terminal_atoms.append(idx)
            if idx in cc_idx:
                neighbors = [x.GetIdx() for x in atom.GetNeighbors()]
                if any([n not in cc_idx for n in neighbors]):
                    last_real_atoms.append(idx)

        logger.info(f"Terminal atoms: {str(list(set(terminal_atoms)))}")
        logger.info(f"Last real atoms: {str(list(set(last_real_atoms)))}")
        
        return (list(set(terminal_atoms)), list(set(last_real_atoms)))

    def _mutate_to_common_core(self, name: str, cc_idx: list, nr_of_steps_for_el: int) -> list:
        """
        Helper function - do not call directly.
        Generates the mutation route to the common fore for mol.
        """
        mol = self.mols[name]
        hydrogens = []
        charge_mutations = []
        lj_mutations = []
        atoms_to_be_mutated = []
        
        # find the atom that connects the common core to the dummy regiom
        terminal_atoms, last_real_atom = self._find_terminal_atom(cc_idx, mol)
        last_real_atom = last_real_atom[0]
        
        # iterate through atoms and select atoms that need to be mutated 
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in cc_idx:
                # 
                if atom.GetSymbol() == 'H' and idx not in terminal_atoms:
                    hydrogens.append(idx)
                atoms_to_be_mutated.append(idx)
                logger.info('Will be decoupled: Idx:{} Element:{}'.format(idx, atom.GetSymbol()))

        if atoms_to_be_mutated:
            ############################################
            ############################################
            # scale all charges of all atoms to zero
            charge_mutations.append(ChargeToZeroMutation(atom_idx=atoms_to_be_mutated,
                                                nr_of_steps=nr_of_steps_for_el, common_core=cc_idx, last_real_atom=last_real_atom))

            ############################################
            ############################################
            # finished with charge mutation
            ############################################
            ############################################
            
            ############################################
            ############################################
            # scale LJ
            ############################################
            ############################################
            # save the last mutation steps
            lj_terminal_mutations = []

            # start with mutation of LJ of hydrogens
            # Only take hydrogens that are not terminal hydrogens
            if hydrogens:
                lj_mutations.append(StericToZeroMutation(hydrogens))
            already_mutated = [] 
            # continue with scaling of heavy atoms LJ
            all_bonds = []
            
            # get all bonds
            for bond in nx.dfs_edges(self.graphs[name]):
                logger.debug(bond)
                all_bonds.append(bond)
            
            for idx1, idx2 in all_bonds:
                # continue if atom is not a hydrogen/already mutated and in the list of to be mutated atoms 
                if idx1 in atoms_to_be_mutated and idx1 not in hydrogens and idx1 not in already_mutated:
                    # is it a terminal atom?
                    if idx1 in terminal_atoms:
                        lj_terminal_mutations.append(StericToDefaultMutation([idx1]))
                    else:
                        lj_mutations.append(StericToZeroMutation([idx1]))
                    already_mutated.append(idx1)
                # continue if atom is not a hydrogen/already mutated and in the list of to be mutated atoms 
                if idx2 in atoms_to_be_mutated and idx2 not in hydrogens and idx2 not in already_mutated:
                    if idx2 in terminal_atoms:
                        lj_terminal_mutations.append(StericToDefaultMutation([idx2]))
                    else:
                        lj_mutations.append(StericToZeroMutation([idx2]))
                    already_mutated.append(idx2)
            
            
            # test that all mutations are included
            # TODO: test that all mutations are covered
            mutations = charge_mutations + lj_mutations + lj_terminal_mutations

            for m in mutations:
                if type(m) == ChargeMutation:
                    logger.debug(f"charge mutation on: {str(m.atom_idx)}")
                elif type(m) == StericMutation:
                    logger.debug(f"steric mutation on: {str(m.atom_idx)}")
                else:
                    logger.debug(f"mutation on: {str(m.atom_idx)}")
        else:
            logger.info("No atoms will be decoupled.")
            mutations = []
        return mutations


class CommonCoreTransformation(object):

    def __init__(self, 
                cc1_indicies: list, 
                cc2_indicies: list, 
                ligand1_psf: pm.charmm.CharmmPsfFile, 
                ligand2_psf: pm.charmm.CharmmPsfFile, 
                nr_of_steps: int, 
                tlc_cc1: str, 
                tlc_cc2: str,
                terminal_atom_idx_cc1: int,
                terminal_atom_idx_cc2: int,
                charge_compensated_ligand2_psf: pm.charmm.CharmmPsfFile,
                charge_mutation:bool,
                bonded_terms_mutation:bool
                ):
        """
        Scale the bonded parameters inside the common core.
        Parameters
        ----------
        cc1_indicies : list
            indices of cc1
        cc2_indicies : list
            indices of cc2 (in the same order as cc1)
        ligand1_psf : pm.charmm.CharmmPsfFile (copy of only ligand)
        ligand2_psf : pm.charmm.CharmmPsfFile (copy of only ligand)
            the target psf that is used to generate the new bonded parmaeters
        nr_of_steps : int
        tlc_cc1 : str
            three letter code of ligand in cc1
        tlc_cc2 : str
            three letter code of ligand in cc2
        """
        self.cc1_indicies = cc1_indicies
        self.cc2_indicies = cc2_indicies
        self.ligand2_psf = ligand2_psf
        self.ligand1_psf = ligand1_psf
        self.nr_of_steps = nr_of_steps
        assert(self.nr_of_steps >= 2)
        self.tlc_cc1 = tlc_cc1
        self.tlc_cc2 = tlc_cc2
        self.terminal_atom_idx_cc1 = terminal_atom_idx_cc1
        self.terminal_atom_idx_cc2 = terminal_atom_idx_cc2
        self.atom_names_mapping, self.terminal_names_mapping = self._get_atom_mapping()
        self.atom_names_mapping_for_bonded_terms = {**self.atom_names_mapping, **self.terminal_names_mapping}
        self.charge_mutation = charge_mutation
        self.bonded_terms_mutation = bonded_terms_mutation
        self.charge_compensated_ligand2_psf = charge_compensated_ligand2_psf
        
        logger.info(f'Bonded terms mutation: {bonded_terms_mutation}')
        logger.info(f'Charge mutation: {charge_mutation}')


    def _get_atom_mapping(self):
        """
        _get_atom_mapping -- match the atom names of the common cores

        Returns
        -------
        [dict]
            matched common core atom names
        """
        # match atomes in common cores
        match_atom_names_cc1_to_cc2 = {}
        for cc1_idx, cc2_idx in zip(self.cc1_indicies, self.cc2_indicies):
            ligand1_atom = self.ligand1_psf[cc1_idx]
            ligand2_atom = self.ligand2_psf[cc2_idx]
            match_atom_names_cc1_to_cc2[ligand1_atom.name] = ligand2_atom.name

        # match terminal atoms
        match_terminal_atoms_cc1_to_cc2 = {self.ligand1_psf[self.terminal_atom_idx_cc1].name : self.ligand2_psf[self.terminal_atom_idx_cc2].name}
        
        return match_atom_names_cc1_to_cc2, match_terminal_atoms_cc1_to_cc2

    def _mutate_charges(self, psf: pm.charmm.CharmmPsfFile, tlc: str, scale: float):

        # common core of psf 1 is transformed to psf 2
        for ligand1_atom in psf.view[f":{tlc}"]:
            if ligand1_atom.name not in self.atom_names_mapping:
                continue
            found = False
            
            # compare to charge compenstated psf 2
            for ligand2_atom in self.charge_compensated_ligand2_psf:
                if self.atom_names_mapping[ligand1_atom.name] == ligand2_atom.name:
                    found = True
                    # are the atoms different?
                    logger.debug(f"Modifying atom: {ligand1_atom}")
                    logger.debug(f"Template atom: {ligand2_atom}")

                    # scale epsilon
                    logger.debug(f"Real charge: {ligand1_atom.charge}")
                    modified_charge = (1.0 - scale) * ligand1_atom.initial_charge + scale * ligand2_atom.charge
                    logger.debug(f"New epsilon: {modified_charge}")
                    ligand1_atom.charge = modified_charge

            if not found:
                raise RuntimeError('No corresponding atom in cc2 found')


    def _mutate_atoms(self, psf: pm.charmm.CharmmPsfFile, tlc: str, scale: float):
        
        """
        mutate atom types. 

        Raises
        ------
        RuntimeError
            if common core atoms can not be matched
        """
        # what will be changed
        mod_type = namedtuple('Atom', 'epsilon, rmin')
        logger.debug('#######################')
        logger.debug('mutate_atoms')

        # iterate through the atoms of the ligand of system1
        for ligand1_atom in psf.view[f":{tlc}"]:
            # continue if not in atom_names_mapping
            if ligand1_atom.name not in self.atom_names_mapping:
                continue

            found = False
            #iterate through the atoms the ligand of system2
            for ligand2_atom in self.ligand2_psf:
                # is there a match up?
                if self.atom_names_mapping[ligand1_atom.name] == ligand2_atom.name:
                    found = True
                    # are the atoms different?
                    if ligand1_atom.type != ligand2_atom.type:
                        self._modify_type(ligand1_atom, psf)
                        logger.debug(f"Modifying atom: {ligand1_atom}")
                        logger.debug(f"Template atom: {ligand2_atom}")

                        # scale epsilon
                        logger.debug(f"Real epsilon: {ligand1_atom.epsilon}")
                        modified_epsilon = (1.0 - scale) * ligand1_atom.epsilon + scale * ligand2_atom.epsilon
                        logger.debug(f"New epsilon: {modified_epsilon}")

                        # scale rmin
                        logger.debug(f"Real rmin: {ligand1_atom.rmin}")
                        modified_rmin = (1.0 - scale) * ligand1_atom.rmin + scale * ligand2_atom.rmin
                        logger.debug(f"New rmin: {modified_rmin}")

                        ligand1_atom.mod_type = mod_type(modified_epsilon, modified_rmin)

            if not found:
                raise RuntimeError('No corresponding atom in cc2 found')

    def _mutate_bonds(self, psf: pm.charmm.CharmmPsfFile, tlc: str, scale: float):

        logger.debug('#######################')
        logger.debug('mutate_bonds')

        mod_type = namedtuple('Bond', 'k, req')
        for ligand1_bond in psf.view[f":{tlc}"].bonds:

            ligand1_atom1_name = ligand1_bond.atom1.name
            ligand1_atom2_name = ligand1_bond.atom2.name
            # all atoms of the bond must be in cc
            # everything outside the cc are bonded terms between dummies or
            # between real atoms and dummies and we can ignore them for now
            if not all(elem in self.atom_names_mapping_for_bonded_terms for elem in [ligand1_atom1_name, ligand1_atom2_name]):
                continue

            found = False
            for ligand2_bond in self.ligand2_psf.bonds:
                ligand2_atom1_name = ligand2_bond.atom1.name
                ligand2_atom2_name = ligand2_bond.atom2.name
                # all atoms of the bond must be in cc
                if not all(elem in self.atom_names_mapping_for_bonded_terms.values() for elem in [ligand2_atom1_name, ligand2_atom2_name]):
                    continue

                # match the two bonds
                if sorted([self.atom_names_mapping_for_bonded_terms[e] for e in [ligand1_atom1_name, ligand1_atom2_name]]) == sorted([ligand2_atom1_name, ligand2_atom2_name]):
                    found = True
                    # are the bonds different?
                    if sorted([ligand1_bond.atom1.type, ligand1_bond.atom2.type]) == sorted([ligand2_bond.atom1.type, ligand2_bond.atom2.type]):
                        continue
                    logger.debug(f"Modifying bond: {ligand1_bond}")

                    logger.debug(f"Template bond: {ligand2_bond}")
                    logger.debug(f'Original value for k: {ligand1_bond.type.k}')
                    logger.debug(f"Target k: {ligand2_bond.type.k}")
                    new_k = ((1.0 - scale) * ligand1_bond.type.k) + (scale * ligand2_bond.type.k)
                    logger.debug(new_k)

                    modified_k = new_k

                    logger.debug(f"New k: {modified_k}")

                    logger.debug(f"Old req: {ligand1_bond.type.req}")
                    modified_req = ((1.0 - scale) * ligand1_bond.type.req) + (scale * ligand2_bond.type.req)
                    logger.debug(f"Modified bond: {ligand1_bond}")

                    ligand1_bond.mod_type = mod_type(modified_k, modified_req)
                    logger.debug(ligand1_bond.mod_type)

            if not found:
                logger.critical(ligand1_bond)
                raise RuntimeError('No corresponding bond in cc2 found: {}'.format(ligand1_bond))

    def _mutate_angles(self, psf: pm.charmm.CharmmPsfFile, tlc: str, scale: float):

        mod_type = namedtuple('Angle', 'k, theteq')
        for cc1_angle in psf.view[f":{tlc}"].angles:
            ligand1_atom1_name = cc1_angle.atom1.name
            ligand1_atom2_name = cc1_angle.atom2.name
            cc1_a3 = cc1_angle.atom3.name

            # only angles in cc
            if not all(elem in self.atom_names_mapping_for_bonded_terms for elem in [ligand1_atom1_name, ligand1_atom2_name, cc1_a3]):
                continue

            found = False
            for cc2_angle in self.ligand2_psf.angles:
                ligand2_atom1_name = cc2_angle.atom1.name
                ligand2_atom2_name = cc2_angle.atom2.name
                cc2_a3 = cc2_angle.atom3.name
                # only angles in cc
                if not all(elem in self.atom_names_mapping_for_bonded_terms.values() for elem in [ligand2_atom1_name, ligand2_atom2_name, cc2_a3]):
                    continue

                if sorted([self.atom_names_mapping_for_bonded_terms[e] for e in [ligand1_atom1_name, ligand1_atom2_name, cc1_a3]]) == sorted([ligand2_atom1_name, ligand2_atom2_name, cc2_a3]):
                    found = True
                    if sorted([cc1_angle.atom1.type, cc1_angle.atom2.type, cc1_angle.atom3.type]) == \
                            sorted([cc2_angle.atom1.type, cc2_angle.atom2.type, cc2_angle.atom3.type]):
                        continue

                    logger.debug(f"Modifying angle: {cc1_angle}")
                    logger.debug(f"Template bond: {cc2_angle}")
                    logger.debug('Scaling k and theteq')

                    logger.debug(f"Old k: {cc1_angle.type.k}")
                    modified_k = (1.0 - scale) * cc1_angle.type.k + scale * cc2_angle.type.k
                    logger.debug(f"New k: {modified_k}")

                    logger.debug(f"Old k: {cc1_angle.type.theteq}")
                    modified_theteq = (1.0 - scale) * cc1_angle.type.theteq + scale * cc2_angle.type.theteq
                    logging.debug(f"New k: {modified_theteq}")

                    cc1_angle.mod_type = mod_type(modified_k, modified_theteq)

            if not found:
                logger.critical(cc1_angle)
                raise RuntimeError('No corresponding angle in cc2 found')

    def _mutate_torsions(self, psf: pm.charmm.CharmmPsfFile, tlc: str, scale: float):

        mod_type = namedtuple('Torsion', 'phi_k, per, phase, scee, scnb')

        # get all torsions present in initial topology
        for cc1_torsion in psf.view[f":{tlc}"].dihedrals:
            ligand1_atom1_name = cc1_torsion.atom1.name
            ligand1_atom2_name = cc1_torsion.atom2.name
            cc1_a3 = cc1_torsion.atom3.name
            cc1_a4 = cc1_torsion.atom4.name
            # all atoms must be in the cc
            if not all(elem in self.atom_names_mapping_for_bonded_terms for elem in [ligand1_atom1_name, ligand1_atom2_name, cc1_a3, cc1_a4]):
                continue

            # get corresponding torsion types in the new topology
            for cc2_torsion in self.ligand2_psf.dihedrals:
                ligand2_atom1_name = cc2_torsion.atom1.name
                ligand2_atom2_name = cc2_torsion.atom2.name
                cc2_a3 = cc2_torsion.atom3.name
                cc2_a4 = cc2_torsion.atom4.name
                # only torsion in cc
                if not all(elem in self.atom_names_mapping_for_bonded_terms.values() for elem in [ligand2_atom1_name, ligand2_atom2_name, cc2_a3, cc2_a4]):
                    continue

                if sorted([self.atom_names_mapping_for_bonded_terms[e] for e in [ligand1_atom1_name, ligand1_atom2_name, cc1_a3, cc1_a4]]) == sorted([ligand2_atom1_name, ligand2_atom2_name, cc2_a3, cc2_a4]):
                    found = True
                    if sorted([cc1_torsion.atom1.type, cc1_torsion.atom2.type, cc1_torsion.atom3.type, cc1_torsion.atom3.type]) == \
                            sorted([cc2_torsion.atom1.type, cc2_torsion.atom2.type, cc2_torsion.atom3.type, cc2_torsion.atom4.type]):
                        continue

                    mod_types = []
                    if scale <= 0.5:
                        # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2
                        for torsion_t in cc1_torsion.type:
                            modified_phi_k = torsion_t.phi_k * max(((1.0 - scale * 2)), 0.0)
                            mod_types.append(mod_type(modified_phi_k, torsion_t.per, torsion_t.phase,
                                                      torsion_t.scee, torsion_t.scnb))
                    else:
                        # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2
                        for torsion_t in cc2_torsion.type:
                            modified_phi_k = torsion_t.phi_k * max((scale - 0.5) * 2, 0.0)
                            mod_types.append(mod_type(modified_phi_k, torsion_t.per, torsion_t.phase,
                                                      torsion_t.scee, torsion_t.scnb))

                    cc1_torsion.mod_type = mod_types
            if not found:
                logger.critical(cc1_torsion)
                raise RuntimeError('No corresponding torsion in cc2 found')

    def mutate(self, psf: pm.charmm.CharmmPsfFile, tlc: str, current_step: int, verbose:int=0):
        """
        Mutates the bonded parameters of cc1 to cc2.
        Parameters
        ----------
        psf : pm.charmm.CharmmPsfFile
            psf that gets mutated
        tlc : str
        current_step : int
            the current step in the mutation protocoll
        only_charge : bool
            only charge is scaled from cc1 to cc2
        """

        assert(type(psf) == pm.charmm.CharmmPsfFile)
        scale = current_step / (self.nr_of_steps)
        if self.charge_mutation:
            logger.info(f" -- Charge parameters from cc1 are transformed to cc2.")
            logger.info(f"Scaling factor:{scale}")
            # scale charge
            self._mutate_charges(psf, tlc, scale)
        elif self.bonded_terms_mutation:
            logger.info(f" -- Atom/Bond/Angle/Torsion parameters from cc1 are transformed to cc2.")
            logger.info(f"Scaling factor:{scale}")
            # scale atoms
            self._mutate_atoms(psf, tlc, scale)
            # scale bonds
            self._mutate_bonds(psf, tlc, scale)
            # scale angles
            self._mutate_angles(psf, tlc, scale)
            # scale torsions
            self._mutate_torsions(psf, tlc, scale)
        else:
            logger.critical('Nothing to do. Is there someting wrong?')

    def _modify_type(self, atom: pm.Atom, psf: pm.charmm.CharmmPsfFile):

        if (hasattr(atom, 'initial_type')):
            # only change parameters
            pass
        else:
            logger.info(f"Setting RRR atomtype for atom: {atom}.")
            atom.type = f"RRR{psf.number_of_dummys}"
            atom.initial_type = atom.type
            psf.number_of_dummys += 1


class BaseMutation(object):

    def __init__(self, atom_idx: list, nr_of_steps: int):
        assert(type(atom_idx) == list)
        self.atom_idx = atom_idx
        self.nr_of_steps = nr_of_steps


class ChargeMutation(BaseMutation):

    def __init__(self, atom_idx: list, nr_of_steps: int, common_core: list):
        if (nr_of_steps <= 2):
            logger.warning('Rapid charge change detected. Proceed with caution.')
        super().__init__(atom_idx, nr_of_steps)
        self.common_core = common_core

    def _scale_charge(self, atom, new_charge):
        """
        Scale atom with charge and compensate the charge change.
        """
        atom.charge = new_charge

    def _compensate_charge(self, psf, tlc:str, old_total_charge:int, last_real_atom:int):
        """Compensate charge change .

        Args:
        """        

        new_charge = round(sum([a.charge for a in psf[f":{tlc.upper()}"].atoms]), 8)
        logger.info('##############')
        logger.info(f"Charge to compensate: {old_total_charge-new_charge}")
        logger.info(f"Adding to atom idx: {psf[last_real_atom]}")
        logger.info('##############')
        
        psf[last_real_atom].charge = psf[last_real_atom].charge+ (old_total_charge-new_charge)
        new_charge = round(sum([a.charge for a in psf[f":{tlc.upper()}"].atoms]),8)

        if not (np.isclose(new_charge, old_total_charge, rtol=1e-4)):
            raise RuntimeError(f'Charge compensation failed. Introducing non integer total charge: {new_charge}.')

        

class StericMutation(BaseMutation):

    def __init__(self, atom_idx: list, nr_of_steps: int):
        super().__init__(atom_idx, nr_of_steps)

    def _scale_epsilon(self, atom, multiplicator):
        logger.debug(atom)
        logger.debug(atom.initial_epsilon)
        atom.epsilon = atom.initial_epsilon * multiplicator

    def _scale_rmin(self, atom, multiplicator):
        logger.debug(atom)
        logger.debug(atom.initial_rmin)
        atom.rmin = atom.initial_rmin * multiplicator

    def _modify_type(self, atom, psf):

        if (hasattr(atom, 'initial_type')):
            # only change parameters
            pass
        else:
            atom.initial_type = atom.type
            atom.type = f"DDD{psf.number_of_dummys}"
            psf.number_of_dummys += 1


class ChargeToZeroMutation(ChargeMutation):

    def __init__(self, atom_idx: list, nr_of_steps: int, common_core: list, last_real_atom:int):
        """
        Scale the electrostatics of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_idx : list
            atoms for which charges are scaled to zero
        mr_of_steps : int
            determines the scaling factor : multiplicator = 1 - (current_step / (nr_of_steps))
        common_core : list
        """
        super().__init__(atom_idx, nr_of_steps, common_core)
        self.last_real_atom = last_real_atom

    def __str__(self):
        return "charges to zero mutation"

    def mutate(self, psf: pm.charmm.CharmmPsfFile, tlc: str, current_step: int):

        """ Performs the mutation """

        old_total_charge = int(round(sum([a.charge for a in psf[f":{tlc.upper()}"].atoms])))
        offset = min([a.idx for a in psf.view[f":{tlc.upper()}"].atoms])
        multiplicator = 1 - (current_step / (self.nr_of_steps))

        logger.info(f" -- Charge to zero mutation.")
        logger.info(f"Scaling factor: {multiplicator}")

        # scale the charge of all atoms
        for idx in self.atom_idx:
            odx = idx + offset
            atom = psf[odx]
            logger.info(f"Scale charge on {atom}")
            logger.info(f"Old charge: {atom.initial_charge}")
            new_charge = float(np.round(atom.initial_charge * multiplicator, 4))
            logger.info(f"New charge: {new_charge}")
            atom.charge = new_charge
            
        if multiplicator != 1:
            # compensate for the total change in charge the terminal atom
            self._compensate_charge(psf, tlc, old_total_charge, self.last_real_atom+ offset)
  
  
  
class StericToDefaultMutation(StericMutation):

    def __init__(self, atom_idx: list):
        """
        Set the steric terms of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_list : list
            atoms for which charges are scaled to zero
        """
        super().__init__(atom_idx, 1)

    def __str__(self):
        return "Steric to zero mutation"

    def __unicode__(self):
        return u"Steric to zero mutation"

    def _modify_type(self, atom, psf):

        if (hasattr(atom, 'initial_type')):
            # only change parameters
            pass
        else:
            atom.initial_type = atom.type
            atom.type = f"DXX0"
            psf.number_of_dummys += 1


    def mutate(self, psf, tlc: str, current_step: int):

        """ Performs the actual mutation """

        logger.info(f" -- Steric interactions to default mutation.")
        logger.info(f"Acting on atoms: {self.atom_idx}")
        offset = min([a.idx for a in psf.view[f":{tlc.upper()}"].atoms])

        for i in self.atom_idx:
            atom = psf[i + offset]
            self._modify_type(atom, psf)
            atom.rmin = 1.5
            atom.epsilon = -0.15
            
class StericToZeroMutation(StericMutation):

    def __init__(self, atom_idx: list):
        """
        Set the steric terms of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_list : list
            atoms for which charges are scaled to zero
        """
        super().__init__(atom_idx, 1)

    def __str__(self):
        return "Steric to zero mutation"

    def __unicode__(self):
        return u"Steric to zero mutation"

    def mutate(self, psf, tlc: str, current_step: int):

        """ Performs the actual mutation """

        logger.info(f" -- Steric interactions to zero mutation.")
        logger.info(f"Acting on atoms: {self.atom_idx}")
        offset = min([a.idx for a in psf.view[f":{tlc.upper()}"].atoms])

        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 0
            self._modify_type(atom, psf)
            self._scale_epsilon(atom, multiplicator)
            self._scale_rmin(atom, multiplicator)

# TODO: we still don't have any implementation for improper torsions!
# TODO: RemoveImproperMutation()

        # # scale improper
        # for torsion in psf.impropers:
        #     # test if one of the torsions that needs to change
        #     if not hasattr(torsion, 'improper_present_in'):
        #         continue

        #     # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2
        #     if torsion.improper_present_in == 'cc1':
        #         # scale the initial torsions down
        #         torsion.type.psi_k = torsion.initial_psi_k * max(((1.0 - scale * 2)), 0.0)
        #     # torsion present at cc1 needs to be turned fully off starting from self.nr_of_steps/2
        #     elif torsion.improper_present_in == 'cc2':
        #         # scale the new torsinos up
        #         torsion.type.psi_k = torsion.initial_psi_k * max((scale -0.5) * 2, 0.0)
