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
from collections import namedtuple
import numpy as np

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
        for idx in self._substructure_match[name]: # BEWARE: the ordering is important - don't cast set!
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

        drawer = rdMolDraw2D.MolDraw2DSVG(800,800)
        drawer.SetFontSize(0.3)

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
        t = self._transform_common_core(nr_of_steps_for_bonded_parameters, nr_of_steps_for_el)

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

    def _transform_common_core(self, nr_of_steps_for_bonded_parameters:int, nr_of_steps_for_el:int)->list:
        """
        Common Core 1 is transformed to Common core 2. Bonded parameters and charges are adjusted. 
        """
        
        transformations = []

        # test if bonded mutations are necessary
        bonded_mutations_necessary = False
        for cc1, cc2 in zip(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2()):
            # did atom type change? if not don't add BondedMutations
            atom1 = self.psfs['m1'][cc1]
            atom2 = self.psfs['m2'][cc2] 
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
            t = BondedParameterMutation(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2(), self.psfs['m1'], self.psfs['m2'], nr_of_steps_for_bonded_parameters, self.s1_tlc, self.s2_tlc)
            transformations.append(t)
        
        # TODO: Add Charge transformation
        logger.info('Charges at commen core need to be transformed')
        t = TransformChargesToTargetCharge(self.get_common_core_idx_mol1(), self.get_common_core_idx_mol2(), self.psfs['m1'], self.psfs['m2'], nr_of_steps_for_el, self.s1_tlc, self.s2_tlc)
        transformations.append(t)

        
        return transformations


    def _mutate_to_common_core(self, name:str, cc_idx:list, nr_of_steps_for_el:int)->list:
        """
        Helper function - do not call directly.
        Generates the mutation route to the common fore for mol2.
        """
        mol = self.mols[name]
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
        mutations.append(ELtoZeroMutation(atom_idx=atoms_to_be_mutated, nr_of_steps=nr_of_steps_for_el, common_core=cc_idx ))
        
        # scale LJ
        # start with mutation of LJ of hydrogens
        mutations.append(LJtoZeroMutation(hydrogens))
        # continue with scaling of heavy atoms LJ
        for idx in atoms_to_be_mutated:
            if idx not in hydrogens:
                mutations.append(LJtoZeroMutation([idx]))

 
        return mutations


class BondedParameterMutation(object):

    def __init__(self, cc1_idx:list, cc2_idx:list, cc1_psf:pm.charmm.CharmmPsfFile, cc2_psf:pm.charmm.CharmmPsfFile, nr_of_steps:int, tlc_cc1:str, tlc_cc2:str):
        """
        Scale the bonded parameters inside the common core.
        Parameters
        ----------
        cc1_idx : list
            indices of cc1
        cc2_idx : list
            indices of cc2 (in the same order as cc1)
        cc1_psf : pm.charmm.CharmmPsfFile (copy of only ligand)
        cc2_psf : pm.charmm.CharmmPsfFile (copy of only ligand)
            the target psf that is used to generate the new bonded parmaeters
        nr_of_steps : int
        tlc_cc1 : str
            three letter code of ligand in cc1
        tlc_cc2 : str
            three letter code of ligand in cc2
        """
        self.cc1_idx = cc1_idx
        self.cc2_idx = cc2_idx
        self.cc2_psf = cc2_psf 
        self.cc1_psf = cc1_psf 
        self.nr_of_steps = nr_of_steps
        assert(self.nr_of_steps >= 2)
        self.tlc_cc1 = tlc_cc1
        self.tlc_cc2 = tlc_cc2
        self.atom_names_mapping = self._get_atom_mapping() 

    def _get_atom_mapping(self):
        match_atom_names_cc1_to_cc2 = {}
        for cc1, cc2 in zip(self.cc1_idx, self.cc2_idx):
            cc1_a = self.cc1_psf[cc1]
            cc2_a = self.cc2_psf[cc2]
            match_atom_names_cc1_to_cc2[cc1_a.name] = cc2_a.name
        
        return match_atom_names_cc1_to_cc2

    def _mutate_atoms(self, psf:pm.charmm.CharmmPsfFile, tlc:str, scale:float):
        mod_type = namedtuple('Atom', 'epsilon, rmin')
        for cc1_atom in psf.view[f":{tlc}"]:
            if cc1_atom.name not in self.atom_names_mapping:
                continue           
            
            found = False
            for cc2_atom in self.cc2_psf:
                if self.atom_names_mapping[cc1_atom.name] == cc2_atom.name:
                    found = True
                    # are the atoms different?
                    if cc1_atom.type != cc2_atom.type:
                        self._modify_type(cc1_atom, psf)
                        logging.info(f"Modifying atom: {cc1_atom}")            
                        logging.info(f"Template atom: {cc2_atom}")
          
                        # scale epsilon
                        logging.info(f"Real epsilon: {cc1_atom.epsilon}")
                        modified_epsilon = (1.0 - scale) * cc1_atom.epsilon + scale * cc2_atom.epsilon
                        logging.info(f"New epsilon: {modified_epsilon}")
                        
                        # scale rmin
                        logging.info(f"Real rmin: {cc1_atom.rmin}")
                        modified_rmin = (1.0 - scale) * cc1_atom.rmin + scale * cc2_atom.rmin
                        logging.info(f"New rmin: {modified_rmin}")

                        cc1_atom.mod_type = mod_type(modified_epsilon, modified_rmin)
            if not found:
                raise RuntimeError('No corresponding atom in cc2 found')
    
    
    def _mutate_bonds(self, psf:pm.charmm.CharmmPsfFile, tlc:str, scale:float):

        mod_type = namedtuple('Bond', 'k, req')
        print(self.atom_names_mapping)
        for cc1_bond in psf.view[f":{tlc}"].bonds:

            cc1_a1 = cc1_bond.atom1.name
            cc1_a2 = cc1_bond.atom2.name
            # all atoms of the bond must be in cc
            # everything outside the cc are bonded terms between dummies or 
            # between real atoms and dummies and we can ignore them for now
            if not all(elem in self.atom_names_mapping for elem in [cc1_a1, cc1_a2]):
                continue

            print(cc1_bond)

            found = False
            for cc2_bond in self.cc2_psf.bonds:
                cc2_a1 = cc2_bond.atom1.name
                cc2_a2 = cc2_bond.atom2.name
                # all atoms of the bond must be in cc
                if not all(elem in self.atom_names_mapping.values() for elem in [cc2_a1, cc2_a2]):
                    continue
                
                print(cc2_bond)
                # match the two bonds
                if sorted([self.atom_names_mapping[e] for e in [cc1_a1, cc1_a2]]) == sorted([cc2_a1, cc2_a2]):
                    found = True
                    print(cc2_bond)
                    # are the bonds different?
                    if sorted([cc1_bond.atom1.type, cc1_bond.atom2.type]) == sorted([cc2_bond.atom1.type, cc2_bond.atom2.type]):
                        continue
                    logging.info('##############################')
                    logging.info(scale)
                    logging.info(f"Modifying bond: {cc1_bond}")

                    logging.info(f"Template bond: {cc2_bond}")
                    logging.info('Original value for k: {}'.format(cc1_bond.type.k))
                    logging.info(f"Target k: {cc2_bond.type.k}")
                    new_k = ((1.0 - scale) * cc1_bond.type.k)   + (scale * cc2_bond.type.k)
                    logging.info(new_k)
                    
                    modified_k = new_k
                    
                    logging.info(f"New k: {modified_k}")

                    logging.info(f"Old req: {cc1_bond.type.req}")
                    modified_req = ((1.0 - scale) * cc1_bond.type.req) + (scale * cc2_bond.type.req)
                    logging.info(f"Modified bond: {cc1_bond}")

                    cc1_bond.mod_type = mod_type(modified_k, modified_req)
                    logger.info(cc1_bond.mod_type)
            
            if not found:
                print(cc1_bond)
                raise RuntimeError('No corresponding bond in cc2 found: {}'.format(cc1_bond))


    def _mutate_angles(self, psf:pm.charmm.CharmmPsfFile, tlc:str, scale:float):

        mod_type = namedtuple('Angle', 'k, theteq')
        for cc1_angle in psf.view[f":{tlc}"].angles:
            cc1_a1 = cc1_angle.atom1.name
            cc1_a2 = cc1_angle.atom2.name
            cc1_a3 = cc1_angle.atom3.name
            
            # only angles in cc
            if not all(elem in self.atom_names_mapping for elem in [cc1_a1, cc1_a2, cc1_a3]):
                    continue

            found = False
            for cc2_angle in self.cc2_psf.angles:
                cc2_a1 = cc2_angle.atom1.name
                cc2_a2 = cc2_angle.atom2.name
                cc2_a3 = cc2_angle.atom3.name
                # only angles in cc
                if not all(elem in self.atom_names_mapping.values() for elem in [cc2_a1, cc2_a2, cc2_a3]):
                    continue

                if sorted([self.atom_names_mapping[e] for e in [cc1_a1, cc1_a2, cc1_a3]]) == sorted([cc2_a1, cc2_a2, cc2_a3]):
                    found = True
                    if sorted([cc1_angle.atom1.type, cc1_angle.atom2.type, cc1_angle.atom3.type]) == \
                        sorted([cc2_angle.atom1.type, cc2_angle.atom2.type, cc2_angle.atom3.type]):
                        continue

                    logging.info(f"Modifying angle: {cc1_angle}")
                    logging.info(f"Template bond: {cc2_angle}")
                    logging.info('Scaling k and theteq')

                    logging.info(f"Old k: {cc1_angle.type.k}")
                    modified_k      = (1.0 - scale) * cc1_angle.type.k      + scale * cc2_angle.type.k
                    logging.info(f"New k: {modified_k}")
                    
                    logging.info(f"Old k: {cc1_angle.type.theteq}")
                    modified_theteq = (1.0 - scale) * cc1_angle.type.theteq + scale * cc2_angle.type.theteq
                    logging.info(f"New k: {modified_theteq}")

                    cc1_angle.mod_type = mod_type(modified_k, modified_theteq)
            
            if not found:
                raise RuntimeError('No corresponding angle in cc2 found')

    def _mutate_torsions(self, psf:pm.charmm.CharmmPsfFile, tlc:str, scale:float):

        mod_type = namedtuple('Torsion', 'phi_k, per, phase, scee, scnb')

        # get all torsions present in initial topology
        for cc1_torsion in psf.view[f":{tlc}"].dihedrals:
            cc1_a1 = cc1_torsion.atom1.name
            cc1_a2 = cc1_torsion.atom2.name
            cc1_a3 = cc1_torsion.atom3.name
            cc1_a4 = cc1_torsion.atom4.name
            # all atoms must be in the cc
            if not all(elem in self.atom_names_mapping for elem in [cc1_a1, cc1_a2, cc1_a3, cc1_a4]):
                continue
            
            # get corresponding torsion types in the new topology
            for cc2_torsion in self.cc2_psf.dihedrals:
                cc2_a1 = cc2_torsion.atom1.name
                cc2_a2 = cc2_torsion.atom2.name
                cc2_a3 = cc2_torsion.atom3.name
                cc2_a4 = cc2_torsion.atom4.name
                # only torsion in cc
                if not all(elem in self.atom_names_mapping.values() for elem in [cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
                    continue

                if sorted([self.atom_names_mapping[e] for e in [cc1_a1, cc1_a2, cc1_a3, cc1_a4]]) == sorted([cc2_a1, cc2_a2, cc2_a3, cc2_a4]):
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
                            modified_phi_k = torsion_t.phi_k * max((scale -0.5) * 2, 0.0)
                            mod_types.append(mod_type(modified_phi_k, torsion_t.per, torsion_t.phase, 
                                                torsion_t.scee, torsion_t.scnb))
                
                    cc1_torsion.mod_type = mod_types
            if not found:
                raise RuntimeError('No corresponding torsion in cc2 found')

    def mutate(self, psf:pm.charmm.CharmmPsfFile, tlc:str, current_step:int):
        """
        Mutates the bonded parameters of cc1 to cc2.
        Parameters
        ----------
        psf : pm.charmm.CharmmPsfFile
            psf that gets mutated
        tlc : str
        current_step : int
            the current step in the mutation protocoll
        """

        assert(type(psf) == pm.charmm.CharmmPsfFile)

        scale = current_step/(self.nr_of_steps)
        # scale atoms
        self._mutate_atoms(psf, tlc, scale)
        # scale bonds
        self._mutate_bonds(psf, tlc, scale)
        # scale angles
        self._mutate_angles(psf, tlc, scale)
        # scale torsions
        self._mutate_torsions(psf, tlc, scale)
            

    def _modify_type(self, atom:pm.Atom, psf:pm.charmm.CharmmPsfFile):

        if (hasattr(atom, 'initial_type')):
            # only change parameters
            pass
        else:
            logging.info(f"Setting RRR atomtype for atom: {atom}.")
            atom.type = f"RRR{psf.number_of_dummys}"
            atom.initial_type = atom.type
            psf.number_of_dummys += 1


class BaseMutation(object):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        assert(type(atom_idx) == list)
        self.atom_idx = atom_idx
        self.nr_of_steps = nr_of_steps


class ELMutation(BaseMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int, common_core:list):
        assert(nr_of_steps >= 2)
        super().__init__(atom_idx, nr_of_steps)
        self.common_core = common_core

    def _scale_charge(self, atom, new_charge):
        """
        Scale atom with charge and compensate the charge change.
        """
        diff_charge = atom.charge - new_charge
        atom.charge = new_charge
        return diff_charge
    
    def _compensate_charge(self, psf, offset:int, diff_charge:float, tlc:str):

        nr_of_atoms_to_spread_charge = len(self.common_core)
        print('##############')
        print('Charge to compensate: {}'.format(diff_charge))
        charge_part = diff_charge / nr_of_atoms_to_spread_charge
        print(charge_part)
        print('##############')
        for idx in self.common_core:
            odx = idx + offset
            psf[odx].charge += charge_part

   
class LJMutation(BaseMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int):
        super().__init__(atom_idx, nr_of_steps)


    def _scale_epsilon(self, atom, multiplicator):
        atom.epsilon = atom.initial_epsilon * multiplicator

    def _scale_rmin(self, atom, multiplicator):
        atom.rmin = atom.initial_rmin *  multiplicator

    def _modify_type(self, atom, psf):

        if (hasattr(atom, 'initial_type')):
            # only change parameters
            pass
        else:
            atom.initial_type = atom.type
            atom.type = f"DDD{psf.number_of_dummys}"
            psf.number_of_dummys += 1

class ELtoZeroMutation(ELMutation):

    def __init__(self, atom_idx:list, nr_of_steps:int, common_core:list):
        """
        Scale the electrostatics of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_list : list
            atoms for which charges are scaled to zero
        mr_of_steps : int
            determines the scaling factor : multiplicator = 1 - (current_step / (nr_of_steps))
        common_core : list
        """
        super().__init__(atom_idx, nr_of_steps, common_core)

    def __str__(self):
        return "charges to zero mutation"
    def __unicode__(self):
        return u"charges to zero mutation"


    def mutate(self, psf:pm.charmm.CharmmPsfFile, tlc:str, current_step:int):
        logger.info('Charges to zero mutation')
        
        old_total_charge = round(sum([a.charge for a in psf[f":{tlc.upper()}"].atoms]))
        offset = min([a.idx for a in psf.view[f":{tlc.upper()}"].atoms])
        diff_charge = 0
        multiplicator = 1 - (current_step / (self.nr_of_steps))
        for idx in self.atom_idx:
            odx = idx + offset
            atom = psf[odx]
            new_charge = float(np.round(atom.initial_charge * multiplicator , 4))
            diff_charge += self._scale_charge(atom, new_charge)
        
        if multiplicator != 1:
            # compensate for the total change in charge 
            self._compensate_charge(psf, offset, diff_charge, tlc)
            new_charges = [a.charge for a in psf[f":{tlc.upper()}"].atoms]
            try:
                assert(np.isclose(sum(new_charges), old_total_charge))
            except AssertionError:
                new_charges = (np.array(new_charges) - np.average(np.array(new_charges)))
                for a, new_charge in zip(psf.view[f":{tlc.upper()}"].atoms, new_charges):
                    a.charge = new_charge
                new_charges = [a.charge for a in psf[f":{tlc.upper()}"].atoms]
                assert(np.isclose(sum(new_charges), old_total_charge))
            

class LJtoZeroMutation(LJMutation):

    def __init__(self, atom_idx:list):
        """
        Set the LJ terms of atoms specified in the atom_idx list to zero.
        Parameters
        ----------
        atom_list : list
            atoms for which charges are scaled to zero
        """
        super().__init__(atom_idx, 1)
    
    def __str__(self):
        return "LJ to zero mutation"
    def __unicode__(self):
        return u"LJ to zero mutation"


    def mutate(self, psf, tlc:str, current_step:int):

        logger.info('LJ to zero mutation')
        offset = min([a.idx for a in psf.view[f":{tlc.upper()}"].atoms])

        for i in self.atom_idx:
            atom = psf[i + offset]
            multiplicator = 0
            self._modify_type(atom, psf)
            self._scale_epsilon(atom, multiplicator)
            self._scale_rmin(atom, multiplicator)
            

class TransformChargesToTargetCharge():
    
    def __init__(self, cc1_idx:list, cc2_idx:list, cc1_psf:pm.charmm.CharmmPsfFile, cc2_psf:pm.charmm.CharmmPsfFile, nr_of_steps:int, tlc_cc1:str, tlc_cc2:str):
        """
        Scale the charges inside the common core.
        Parameters
        ----------
        cc1_idx : list
            indices of cc1
        cc2_idx : list
            indices of cc2 (in the same order as cc1)
        cc1_psf : pm.charmm.CharmmPsfFile (copy of only ligand)
        cc2_psf : pm.charmm.CharmmPsfFile (copy of only ligand)
            the target psf that is used to generate the new bonded parmaeters
        nr_of_steps : int
        tlc_cc1 : str
            three letter code of ligand in cc1
        tlc_cc2 : str
            three letter code of ligand in cc2
        """
        self.cc1_idx = cc1_idx
        self.cc2_idx = cc2_idx
        self.cc2_psf = cc2_psf 
        self.cc1_psf = cc1_psf 
        self.tlc_cc1 = tlc_cc1
        self.tlc_cc2 = tlc_cc2
        assert(nr_of_steps >= 2)
        self.nr_of_steps = nr_of_steps
        self.atom_names_mapping = self._get_atom_mapping()
    
    def __str__(self):
        return "Transform charge distribution of common core 1 to common core 2"
    def __unicode__(self):
        return u"Transform charge distribution of common core 1 to common core 2"

    def _get_atom_mapping(self):
        match_atom_names_cc1_to_cc2 = {}
        for cc1, cc2 in zip(self.cc1_idx, self.cc2_idx):
            cc1_a = self.cc1_psf[cc1]
            cc2_a = self.cc2_psf[cc2]
            match_atom_names_cc1_to_cc2[cc1_a.name] = cc2_a.name
        
        return match_atom_names_cc1_to_cc2

    def _compensate_charge(self, psf, diff_charge:float, total_charge):

        nr_of_atoms_to_spread_charge = len(self.cc2_idx)
        print('##############')
        print('Charge to compensate: {}'.format(diff_charge))
        charge_part = diff_charge / nr_of_atoms_to_spread_charge
        print('##############')
        for idx in self.cc2_idx:
            psf[idx].charge += charge_part

        return psf

    def _scale_cc2_charges(self):
        """ set all charges not in cc to zero"""
        
        cc2_scaled_psf_ligand = self.cc2_psf[f":{self.tlc_cc2.upper()}"] 
        new_charge = 0.0
        diff_charge = 0.0
        for atom in cc2_scaled_psf_ligand:
            if atom.idx in self.cc2_idx:
                continue
            diff_charge += atom.charge - new_charge    
            atom.charge = new_charge
        return cc2_scaled_psf_ligand, diff_charge



    def _mutate_charge(self, psf:pm.charmm.CharmmPsfFile, tlc:str, current_step:int):
        """ mutate charges of cc1 to cc2"""
        
        scale =  (current_step / (self.nr_of_steps))
        total_charge = round(sum([a.charge for a in self.cc2_psf[f":{self.tlc_cc2.upper()}"].atoms]))
        cc2_scaled_psf_ligand, diff_charge = self._scale_cc2_charges()
        cc2_psf = self._compensate_charge(cc2_scaled_psf_ligand, diff_charge, total_charge)

        for cc1_atom in psf.view[f":{tlc}"]:
            if cc1_atom.name not in self.atom_names_mapping:
                continue
            
            logging.info(f"Scale charge for atom: {cc1_atom}")
            
            for cc2_atom in cc2_psf:
                if self.atom_names_mapping[cc1_atom.name] == cc2_atom.name:
                    break

            logging.info(f"Template atom: {cc2_atom}")
            # scale charge # NOTE: Charges are scaled directly
            logging.info(f"Old charge: {cc1_atom.charge}")
            
            modified_charge = (1.0 - scale) * cc1_atom.initial_charge + scale * cc2_atom.charge

            cc1_atom.charge = modified_charge
            logging.info(f"New charge: {cc1_atom.charge}")



    def mutate(self, psf:pm.charmm.CharmmPsfFile, tlc:str, current_step:int):
        """
        Mutates the charge parameters of cc1 to cc2.
        Parameters
        ----------
        psf : pm.charmm.CharmmPsfFile
            psf that gets mutated
        tlc : str
        """

        assert(type(psf) == pm.charmm.CharmmPsfFile)

        # mutate charge
        self._mutate_charge(psf, tlc, current_step)



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




