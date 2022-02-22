import logging

import networkx as nx
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS

logger = logging.getLogger(__name__)

"""
General function for finding the common core of two molecules, creating dictionaries, converting rdkit-mol to networkx-graph etc. 
"""


def generate_apply_dicts(mol):
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

    for atom in mol.GetAtoms():
        atom_name = atom.GetSymbol()
        #  atom_name = atom.GetAtomicNum()
        atom_index = atom.GetIdx()
        atom_type = atom.GetSymbol()
        atom_charge = atom.GetFormalCharge()

        atom_idx_to_atom_name[atom_index] = atom_name
        atom_name_to_atom_idx[atom_name] = atom_index
        atom_name_to_atom_type[atom_name] = atom_type
        atom_idx_to_atom_partial_charge[atom_index] = atom_charge

    for atom in mol.GetAtoms():

        atom.SetProp("atom_name", atom_idx_to_atom_name[atom.GetIdx()])
        atom.SetProp(
            "atom_type",
            atom_name_to_atom_type[atom_idx_to_atom_name[atom.GetIdx()]],
        )
        atom.SetProp("atom_index", str(atom.GetIdx()))
        atom.SetProp("atom_charge", str(atom_idx_to_atom_partial_charge[atom.GetIdx()]))
    return mol


# ### diverse functions for converting rdkit-mol object to networkx-graph object


def mol_to_nx(mol):
    """
    function for converting rdkit-mol object to networkx-graph object
    without creating weights
    therefore not used for new mutation algorithms (only for comparison with old dfs-algorithm)
    ----
    returns a nx-graph-object
    """
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            is_aromatic=atom.GetIsAromatic(),
            atom_symbol=atom.GetSymbol(),
        )

    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType()
        )

    return G


def _mol_to_nx_full(mol):
    """
    function for converting rdkit-mol object to networkx-graph object
    without creating weights
    therefore not used for new mutation algorithms (only for comparison with old dfs-algorithm)
    several further attributes are added to get same representation in networkx as for the original Transformato-dfs-mutation-algorithm (solely for comparison)
    ----
    returns a nx-graph-object
    """
    G = nx.Graph()

    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp(
            "molAtomMapNumber2", str(mol.GetAtomWithIdx(idx).GetIdx())
        )

    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
            atomic_id=atom.GetIdx(),
            atom_symbol=atom.GetSymbol(),
            atom_map=atom.GetProp("molAtomMapNumber2"),
            # additionally
            atom_index=atom.GetProp("atom_index"),
            atom_type=atom.GetProp("atom_type"),
            atom_index_type=str(atom.GetProp("atom_type"))
            + ":"
            + str(atom.GetProp("atom_index"))
            #  str(i.GetProp("atom_index")) + ":" + i.GetProp("atom_type")
        )

    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType(),
        )
    return G


def _mol_to_nx_full_weight(
    mol, indiv_atom_weights: bool = False, atom_weights: dict = {"N": 50, "C": 50}
):
    """
    function for converting rdkit-mol object to networkx-graph object
    also weights are added (necessary for new mutation algorithms)
    ----
    returns a nx-graph-object with weight attribute
    """

    G = nx.Graph()

    # to get same representation in networkx as in transformato(copied from mutate.py)
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp(
            "molAtomMapNumber2", str(mol.GetAtomWithIdx(idx).GetIdx())
        )

    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
            atomic_id=atom.GetIdx(),
            atom_symbol=atom.GetSymbol(),
            atom_map=atom.GetProp("molAtomMapNumber2"),
            # additionally
            atom_index=atom.GetProp("atom_index"),
            atom_type=atom.GetProp("atom_type"),
            atom_index_type=str(atom.GetProp("atom_type"))
            + ":"
            + str(atom.GetProp("atom_index"))
            #  str(i.GetProp("atom_index")) + ":" + i.GetProp("atom_type")
        )

    for bond in mol.GetBonds():

        weight = 50
        if indiv_atom_weights == True:
            importantatom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
            importantatom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
            for key in atom_weights:
                if (
                    key == importantatom1.GetSymbol()
                    or key == importantatom2.GetSymbol()
                ):
                    weight = atom_weights[key]

        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType(),
            weight=weight,
        )
    return G


def get_common_core(mol1, mol2):
    """
    get the common core of two molecules (rdkit-mols)
    """

    mols = [mol1, mol2]

    res = rdFMCS.FindMCS(
        mols,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        ringCompare=rdkit.Chem.rdFMCS.RingCompare.StrictRingFusion,
    )
    substructure = Chem.MolFromSmarts(res.smartsString)
    # res=rdFMCS.FindMCS(mols, ringMatchesRingOnly = True, completeRingsOnly = True )

    ccore = Chem.rdmolfiles.MolFragmentToCXSmiles(
        mol1,
        mol1.GetSubstructMatch(substructure),
        kekuleSmiles=True,
        isomericSmiles=False,
    )
    ccoremol = Chem.MolFromSmiles(ccore)

    hit_ats1 = mol1.GetSubstructMatch(substructure)
    hit_bonds1 = []
    for bond in substructure.GetBonds():
        aid1 = hit_ats1[bond.GetBeginAtomIdx()]
        aid2 = hit_ats1[bond.GetEndAtomIdx()]
        hit_bonds1.append(mol1.GetBondBetweenAtoms(aid1, aid2).GetIdx())

    hit_ats2 = mol2.GetSubstructMatch(substructure)
    hit_bonds2 = []
    for bond in substructure.GetBonds():
        aid1 = hit_ats2[bond.GetBeginAtomIdx()]
        aid2 = hit_ats2[bond.GetEndAtomIdx()]
        hit_bonds2.append(mol2.GetBondBetweenAtoms(aid1, aid2).GetIdx())

    mol1coreindex = [hit_ats1]
    mol2coreindex = [hit_ats2]

    return mol1coreindex, mol2coreindex, hit_ats1, hit_ats2


def _find_connected_dummy_regions_mol(mol, ccore, G: nx.Graph):
    """
    find connected dummy regions
    ---
    mol: rdkit-mol of molecule
    ccore: previously computed ccore for molecule
    G: networkx-Graph representation of molecule
    """

    sub = ccore
    G_dummy = G.copy()

    from itertools import chain

    sub = [int(x) for x in chain.from_iterable(ccore)]

    # find all dummy atoms
    list_of_dummy_atoms_idx = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in sub
    ]
    nr_of_dummy_atoms = len(list_of_dummy_atoms_idx) + 1
    list_of_real_atoms_idx = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() in sub
    ]

    logger.info("list_of_dummy_atoms_idx")
    logger.info(list_of_dummy_atoms_idx)

    logger.info("nr_of_dummy_atoms")
    logger.info(nr_of_dummy_atoms)

    logger.info("list_of_real_atoms_idx")
    logger.info(list_of_real_atoms_idx)

    # remove real atoms from graph to obtain multiple connected compounds
    for real_atom_idx in list_of_real_atoms_idx:
        G_dummy.remove_node(real_atom_idx)

    # find these connected compounds
    from networkx.algorithms.components import connected_components

    unique_subgraphs = [
        c for c in sorted(nx.connected_components(G_dummy), key=len, reverse=True)
    ]

    return unique_subgraphs, G_dummy


def _find_terminal_atom(cc_idx: int, mol):
    """
    Find atoms that connect the molecule to the common core.
    Args:
        cc_idx (list): common core index atoms
        mol ([type]): rdkit mol object
    """
    terminal_dummy_atoms = []
    terminal_real_atoms = []

    from itertools import chain

    cc_idx = [int(x) for x in chain.from_iterable(cc_idx)]

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()

        if idx not in cc_idx:
            neighbors = [x.GetIdx() for x in atom.GetNeighbors()]
            if any([n in cc_idx for n in neighbors]):
                terminal_dummy_atoms.append(idx)
        if idx in cc_idx:
            neighbors = [x.GetIdx() for x in atom.GetNeighbors()]
            if any([n not in cc_idx for n in neighbors]):
                terminal_real_atoms.append(idx)

    logger.info(f"Terminal dummy atoms: {str(list(set(terminal_dummy_atoms)))}")
    logger.info(f"Terminal real atoms: {str(list(set(terminal_real_atoms)))}")

    return (list(set(terminal_dummy_atoms)), list(set(terminal_real_atoms)))


def _match_terminal_real_and_dummy_atoms(
    mol, real_atoms_cc: list, dummy_atoms_cc: list
):
    """
    Matches the terminal real and dummy atoms and returns a dict with real atom idx as key and a set of dummy atoms that connect
    to this real atom as a set
    Parameters
    ----------
    mol : [Chem.Mol]
        The mol object with the real and dummy atoms
    real_atoms_cc : list
        list of real atom idx
    dummy_atoms_cc : list
        list of dummy atom idx

    Returns
    -------
    [type]
        [description]
    """

    from collections import defaultdict

    real_atom_match_dummy_atom = defaultdict(set)
    for real_atom_idx in real_atoms_cc:
        real_atom = mol.GetAtomWithIdx(real_atom_idx)
        real_neighbors = [x.GetIdx() for x in real_atom.GetNeighbors()]
        for dummy_atoms_idx in dummy_atoms_cc:
            if dummy_atoms_idx in real_neighbors:
                real_atom_match_dummy_atom[real_atom_idx].add(dummy_atoms_idx)

    return real_atom_match_dummy_atom


def reduce_terminal(match_terminal_atoms: dict, subg, G_dummy: nx.Graph):
    """
    again check for connected components after matching real and dummy atoms
    if atom appears in multiple components, remove one of the connections
    (only necessary for pathological/illegal cases where two atoms of the same component are connected to the common core therefore should not be used in regular workflow)
    """
    coreconnectingnodes = []

    for nodes in match_terminal_atoms.values():

        for node in nodes:
            coreconnectingnodes.append(node)

        for i in range(len(coreconnectingnodes)):

            for j in range(i, len(coreconnectingnodes)):

                if i != j:

                    if nx.has_path(
                        G_dummy, coreconnectingnodes[i], coreconnectingnodes[j]
                    ):

                        nodes.remove(coreconnectingnodes[i])

        coreconnectingnodes = []

    aequinodes = []

    for nodes in match_terminal_atoms.values():
        for node in nodes:
            aequinodes.append(node)

    removenodes = []
    for i in range(len(aequinodes)):
        for j in range(len(aequinodes)):
            if i != j:
                for elem in subg:
                    if aequinodes[i] in elem and aequinodes[j] in elem:
                        if (
                            aequinodes[i] not in removenodes
                            and aequinodes[j] not in removenodes
                        ):
                            removenodes.append(aequinodes[i])

    for removenode in removenodes:

        match_terminal_atoms = {
            k: v for k, v in match_terminal_atoms.items() if removenode not in v
        }

    return match_terminal_atoms
