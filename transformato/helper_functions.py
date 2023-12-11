import logging
import parmed as pm
import networkx as nx
from openmm import unit
import sys

logger = logging.getLogger(__name__)


#### Constants which might be needed sometimes!
temperature = 303.15 * unit.kelvin
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kT = kB * temperature
charmm_gpu = ""
# charmm_gpu = "domdec-gpu"  # uncomment this if you want to use domdec-gpu
# charmm_gpu = "charmm-openmm"  # uncomment this if you want to use charmm-openmm
this = sys.modules[__name__]
# we can explicitly make assignments on it
this.NUM_PROC = 1



#####################################
# config for tets
test_platform_openMM = "GPU"
test_platform_CHARMM = "CPU"
test_platform_override = True
#####################################


def change_platform_to_test_platform(configuration: dict, engine: str):
    if engine == "openMM":
        change_to = test_platform_openMM
    elif engine == "CHARMM":
        change_to = test_platform_CHARMM
    else:
        raise RecursionError()

    if change_to.upper() == "GPU":
        configuration["simulation"]["GPU"] = True
        print("Setting platform to GPU")
    elif change_to.upper() == "CPU":
        configuration["simulation"]["GPU"] = False
        print("Setting platform to CPU")
    else:
        raise RuntimeError("something went wrong")


def calculate_order_of_LJ_mutations_asfe(central_atoms: list, G: nx.Graph) -> list:
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


def change_route_cycles(route, cycledict, degreedict, weightdict, G):
    """
    preliminary mutation list is sorted using cycle and degree dictionary
    currently used in _calculate_order_of_LJ_mutations_new
    ----
    Args:
        route: original mutation route
        cycledict: dict of cycle participation of atoms
        degreedict: dict of  degree of atoms
        weightdict: dict of  weight of atoms
        G: nx-graph of molecule
    ----
    returns reordered array of the mutation route
    """

    for i in range(len(route) - 1):
        routedict = route[i]
        routeweight = weightdict.get(route[i])

        routecycleval = cycledict.get(route[i])
        routedegreeval = degreedict.get(route[i])

        for j in range(i, len(route)):
            if routeweight == weightdict[route[j]]:
                # if nodes have same weight (i.e. distance from root), the node participating in more cycles is removed later

                if routecycleval > cycledict[route[j]] or (
                    routecycleval == cycledict[route[j]]
                    and routedegreeval > degreedict[route[j]]
                ):
                    idx1 = route.index(route[i])
                    idx2 = route.index(route[j])
                    route[idx1], route[idx2] = route[idx2], route[idx1]
                    continue

                # if nodes have same weight (i.e. distance from root) and same cycle participation number, the node which has more neighbours already removed is removed earlier

                if routecycleval == cycledict[route[j]]:
                    edgesi = G.edges(routedict)
                    edgesj = G.edges(route[j])

                    iedgecounter = 0
                    for edge in edgesi:
                        if edge[1] in route[0:i] or edge[0] in route[0:i]:
                            iedgecounter = iedgecounter + 1

                    jedgecounter = 0
                    for edge in edgesj:
                        if edge[1] in route[0:i] or edge[0] in route[0:i]:
                            jedgecounter = jedgecounter + 1

                    if iedgecounter < jedgecounter:
                        idx1 = route.index(route[i])
                        idx2 = route.index(route[j])
                        route[idx1], route[idx2] = route[idx2], route[idx1]

    return route


def cycle_checks_nx(G, use_actual_weight_for_mod=False):
    """
    cycle processing, can be used in _calculate_order_of_LJ_mutations_new_iter and ..._new_iter_change (default is cycle_checks_nx_v2)
    --------
    returns nx-graph-object with updated weights (according to cycle participation of the atom)
    """

    # search cycles using networkx
    cycles = nx.cycle_basis(G)

    from collections import Counter

    cdict = Counter(x for xs in cycles for x in set(xs))

    # modify weighted graph: nodes participating in many cycles get lower weight
    for i in cdict:
        edg = G.edges(i)
        for el in edg:
            if use_actual_weight_for_mod == True:
                G[el[0]][el[1]]["weight"] = G[el[0]][el[1]]["weight"] - cdict[i] * (
                    G[el[0]][el[1]]["weight"] ** (1 / 2)
                )
            else:
                G[el[0]][el[1]]["weight"] = G[el[0]][el[1]]["weight"] - cdict[i] * 5

    return G


def cycle_checks(G):
    """
    cycle processing dictionary and degree dictionary for preferential removal (atoms which neighbours already have been removed are removed earlier), currently used in _calculate_order_of_LJ_mutations_new (via change_route_cycles)
    ----
    returns dictionary containing number of cycle participation (cdict) and dict containing degree of atom (degreedict)
    """

    # search cycles using networkx
    cycles = nx.cycle_basis(G)

    # alternatively, using rdkit
    # ri = mol.GetRingInfo()
    # cyclesrdkit = ri.AtomRings()

    import collections
    from collections import Counter

    cdict = Counter(x for xs in cycles for x in set(xs))
    # cdictrdkit = Counter(x for xs in cyclesrdkit for x in set(xs))

    # add atoms with no cycle participation
    for key in G.nodes:
        if key not in cdict:
            cdict[key] = 0

    degreedict = G.degree()
    degreedict = {node: val for (node, val) in degreedict}

    return cdict, degreedict


def exclude_Hs_from_mutations(connected_dummy_regions: list, G: nx.Graph):
    """
    hydrogens are removed from the networkx-graph-representation and the list of connected dummy regions
    ----
    Args:
        connected_dummy_regions: list of connected dummy regions
        G: nx-graph of molecule
    ----
    returns list of connected dummy regions and networkx-graph without hydrogens
    """

    G_hydrogens = [x for x, y in G.nodes(data=True) if y["atom_type"] == "H"]

    G.remove_nodes_from(G_hydrogens)
    connected_dummy_regions_copy = connected_dummy_regions
    for hydroindex in G_hydrogens:
        for indexregion, region in enumerate(connected_dummy_regions):
            if hydroindex in region:
                connected_dummy_regions_copy[indexregion].remove(hydroindex)

    return connected_dummy_regions_copy, G
