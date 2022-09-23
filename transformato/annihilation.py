import logging
from collections import namedtuple, defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np
import networkx as nx
import parmed as pm
from IPython.core.display import display
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.molSize = (900, 900)  # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

from transformato.system import SystemStructure
from transformato.mutate import DummyRegion, MutationDefinition

logger = logging.getLogger(__name__)


class ProposeMutationRouteASFE(object):
    def __init__(
        self,
        s1: SystemStructure,
    ):
        """
        A class that proposes the mutation route for only one molecule.
        Parameters
        ----------
        mol1: Chem.Mol
        """

        mol1_name: str = "m1"
        self.system: dict = {"system1": s1}
        self.mols: dict = {mol1_name: s1.mol}
        self.graphs: dict = {mol1_name: s1.graph}
        # psfs for reference of only ligand
        self.psfs: dict = {s1.psfs["waterbox"][f":{s1.tlc}"]}
        self.psf1: pm.charmm.CharmmPsfFile = s1.psfs
        self._substructure_match: dict = {mol1_name: []}
        self.removed_indeces: dict = {mol1_name: []}
        self.added_indeces: dict = {mol1_name: []}
        self.s1_tlc = s1.tlc

        self.dummy_region_cc1: DummyRegion

    @staticmethod
    def _calculate_order_of_LJ_mutations(central_atoms: list, G: nx.Graph) -> list:

        ordered_LJ_mutations = []
        root = central_atoms[0]

        final_order = []
        removearray = []
        removeG = nx.Graph()
        G_dummy = G.copy()

        while len(G_dummy.nodes()) > 0:

            G_origweights = G_dummy.copy()

            # dijkstra
            ssource = nx.single_source_dijkstra(G_dummy, source=root, weight="weight")

            # find node with maximum distance - i.e next to be removed
            max_node = max(ssource[0], key=ssource[0].get)

            # the currently most distant node is added to the final order - it is the next to be removed
            final_order.extend([max_node])

            # restore original weights
            G_dummy = G_origweights

            # remove the most distant node (after added to the final_order-array); the next iteration will find the most distant node of the reduced graph
            G_dummy.remove_node(max_node)

            # add to removeG
            removeG.add_node(max_node)
            removearray.append(max_node)

            # final_order contains mutation order for current dummy region
            sortedssource_edges = final_order

        # sortedssource_edges_list already contains the nodes in right (reversed) order (including root)
        nodes = sortedssource_edges
        ordered_LJ_mutations.append(nodes)

        logger.info(
            f"Turning off the LJ potential in the following order: {ordered_LJ_mutations}"
        )

        return ordered_LJ_mutations

    def propose_common_core(self):
        self._get_idx_of_all_atoms("m1")

    def finish_common_core(
        self,
        odered_connected_dummy_regions_cc1: list = [],
    ):

        logger.info("Find connected dummy regions")
        # all atoms should become dummy atoms in the end

        central_atoms = nx.center(self.graphs["m1"])

        # Assure, that the central atom is no hydrogen
        for atom in self.psf1["waterbox"][f":{self.s1_tlc}"].atoms:
            if atom.idx in central_atoms:
                if atom.name.startswith("H") == True:
                    raise RuntimeError(
                        f"One of the central atoms seems to be a hydrogen atom"
                    )

        # calculate the ordering or LJ mutations
        if not odered_connected_dummy_regions_cc1:
            odered_connected_dummy_regions_cc1 = self._calculate_order_of_LJ_mutations(
                central_atoms,
                self.graphs["m1"].copy(),
            )

        self.dummy_region_cc1 = DummyRegion(
            mol_name="m1",
            tlc=self.s1_tlc,
            match_termin_real_and_dummy_atoms=[],
            connected_dummy_regions=odered_connected_dummy_regions_cc1,
            lj_default=[],
        )

    def get_common_core_idx_mol1(self) -> list:
        """
        Returns the common core of mol1.
        """
        return self._get_common_core("m1")

    def _get_common_core(self, name: str) -> list:
        """
        Helper Function - should not be called directly.
        Returns the common core.
        """
        keep_idx = []
        # BEWARE: the ordering is important - don't cast set!
        for idx in self._substructure_match[name] + self.added_indeces[name]:
            if idx not in self.removed_indeces[name]:
                keep_idx.append(idx)
        return keep_idx

    def _get_idx_of_all_atoms(
        self,
        mol1_name: str,
    ):
        """
        Iterates over all atoms of the molecule and saves them as a list
        ----------
        mol1_name: str
        """

        s1 = []
        for atom in self.psf1["waterbox"][f":{self.s1_tlc}"].atoms:
            s1.append(atom.idx)

        self._substructure_match[mol1_name] = list(s1)

    def _return_atom_idx_from_bond_idx(self, mol: Chem.Mol, bond_idx: int):
        return (
            mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx(),
            mol.GetBondWithIdx(bond_idx).GetEndAtomIdx(),
        )

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
                mol.GetAtomWithIdx(idx).SetProp(
                    "molAtomMapNumber", str(mol.GetAtomWithIdx(idx).GetIdx())
                )
            return mol

        mol = mol_with_atom_index(mol)
        AllChem.Compute2DCoords(mol)
        display(mol)

    def generate_mutations_to_common_core_for_mol1(self) -> dict:
        """
        Generates the mutation route to the common fore for mol1.
        ----------
        mutations: list
            list of mutations
        """

        m = self._mutate_to_common_core(
            self.dummy_region_cc1, self.get_common_core_idx_mol1(), mol_name="m1"
        )

        return m

    def _mutate_to_common_core(
        self, dummy_region: DummyRegion, cc_idx: list, mol_name: str
    ) -> dict:
        """
        Helper function - do not call directly.
        Generates the mutation route to the common fore for mol.
        """
        mutations = defaultdict(list)
        # iterate through atoms and select atoms that need to be mutated
        atoms_to_be_mutated = []
        hydrogens = []

        for atom in self.psf1["waterbox"][f":{self.s1_tlc}"].atoms:
            # idx = atom.idx - self.offset
            idx = atom.idx
            if idx in dummy_region.connected_dummy_regions[0]:
                if atom.name.startswith("H") == True:
                    hydrogens.append(idx)
                atoms_to_be_mutated.append(idx)
                logger.info(
                    "Will be decoupled: Idx:{} Element:{}".format(idx, atom.name)
                )

        if atoms_to_be_mutated:
            ############################################
            ############################################
            # charge mutation
            ############################################
            ############################################

            m = MutationDefinition(
                atoms_to_be_mutated=atoms_to_be_mutated,
                common_core=[],
                dummy_region=dummy_region,
                vdw_atom_idx=[],
                steric_mutation_to_default=False,
            )

            mutations["charge"].append(m)

            ############################################
            ############################################
            # LJ mutation
            ############################################
            ############################################

            # start with mutation of LJ of hydrogens
            # Only take hydrogens that are not terminal hydrogens

            if hydrogens:
                m = MutationDefinition(
                    atoms_to_be_mutated=atoms_to_be_mutated,
                    common_core=[],
                    dummy_region=dummy_region,
                    vdw_atom_idx=hydrogens,
                    steric_mutation_to_default=False,
                )

                mutations["hydrogen-lj"].append(m)

            for atom_idx in dummy_region.connected_dummy_regions[0]:
                if atom_idx in hydrogens:
                    # is already mutated in the previouse step
                    continue
                else:
                    # normal lj mutation
                    m = MutationDefinition(
                        atoms_to_be_mutated=atoms_to_be_mutated,
                        common_core=[],
                        dummy_region=dummy_region,
                        vdw_atom_idx=[atom_idx],
                        steric_mutation_to_default=False,
                    )
                    mutations["lj"].append(m)

        else:
            logger.critical("No atoms will be decoupled.")
            mutations = defaultdict()

        return mutations

    def show_common_core_on_mol1(self, show_atom_types: bool = False):
        """
        Shows common core on mol1
        """
        return self._show_common_core(
            self.mols["m1"], self.get_common_core_idx_mol1(), show_atom_types
        )

    def show_common_core_on_mol2(self, show_atom_types: bool = False):
        """
        Shows common core on mol2
        """
        return self._show_common_core(
            self.mols["m2"], self.get_common_core_idx_mol2(), show_atom_types
        )

    def _show_common_core(self, mol, highlight: list, show_atom_type: bool):
        """
        Helper function - do not call directly.
        Show common core.
        """
        # https://rdkit.blogspot.com/2015/02/new-drawing-code.html

        mol = deepcopy(mol)
        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(800, 800)
        drawer.SetFontSize(6)

        opts = drawer.drawOptions()

        if show_atom_type:
            for i in mol.GetAtoms():
                opts.atomLabels[i.GetIdx()] = (
                    str(i.GetProp("atom_index")) + ":" + i.GetProp("atom_type")
                )
        else:
            for i in mol.GetAtoms():
                opts.atomLabels[i.GetIdx()] = (
                    str(i.GetProp("atom_index")) + ":" + i.GetProp("atom_name")
                )

        drawer.DrawMolecule(mol, highlightAtoms=highlight)
        Draw.DrawingOptions.includeAtomNumbers = False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")
        return svg
