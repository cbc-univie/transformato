import logging

import networkx as nx

logger = logging.getLogger(__name__)

"""
4 different functions for calculating mutations (_calculate_order_of_LJ_mutations, _calculate_order_of_LJ_mutations_new, _calculate_order_of_LJ_mutations_iter and _calculate_order_of_LJ_mutations_iter_change):
_calculate_order_of_LJ_mutations: naive dfs (as currently in transformato)
_calculate_order_of_LJ_mutations_new: bfs/djikstra-algorithm applied once for route
_calculate_order_of_LJ_mutations_new_iter: bfs/djikstra-algorithm applied iteratively, i.e. after each removal of an atom 
_calculate_order_of_LJ_mutations_new_iter_change: works iteratively, i.e. after each removal of an atom, algorithm is chosen depending on current state
"""


def _calculate_order_of_LJ_mutations(
    connected_dummy_regions: list, match_terminal_atoms: dict, G: nx.Graph
) -> list:
    """
    dfs mutation algorithm (as currently in transformato)
    """

    ordered_LJ_mutations = []
    for real_atom in match_terminal_atoms:
        for dummy_atom in match_terminal_atoms[real_atom]:
            for connected_dummy_region in connected_dummy_regions:
                # stop at connected dummy region with specific dummy_atom in it
                if dummy_atom not in connected_dummy_region:
                    continue

                G_dummy = G.copy()
                # delete all nodes not in dummy region
                remove_nodes = [
                    node for node in G.nodes() if node not in connected_dummy_region
                ]
                for remove_node in remove_nodes:
                    G_dummy.remove_node(remove_node)

                # root is the dummy atom that connects the real region with the dummy region
                root = dummy_atom

                edges = list(nx.dfs_edges(G_dummy, source=root))

                nodes = [root] + [v for u, v in edges]

                nodes.reverse()  # NOTE: reverse the mutation
                ordered_LJ_mutations.append(nodes)

    return ordered_LJ_mutations


def _calculate_order_of_LJ_mutations_new(
    connected_dummy_regions: list,
    match_terminal_atoms: dict,
    G: nx.Graph,
    cyclecheck=True,
    ordercycles=True,
) -> list:
    """
    bfs/djikstra-algorithm applied once for route (without iterations)
    """

    ordered_LJ_mutations = []

    for real_atom in match_terminal_atoms:
        for dummy_atom in match_terminal_atoms[real_atom]:
            for connected_dummy_region in connected_dummy_regions:
                # stop at connected dummy region with specific dummy_atom in it
                if dummy_atom not in connected_dummy_region:
                    continue

                G_dummy = G.copy()
                # delete all nodes not in dummy region
                remove_nodes = [
                    node for node in G.nodes() if node not in connected_dummy_region
                ]
                for remove_node in remove_nodes:
                    G_dummy.remove_node(remove_node)

                # root is the dummy atom that connects the real region with the dummy region
                root = dummy_atom

                # process cycles
                if cyclecheck == True and ordercycles == False:
                    G_dummy = cycle_checks_nx(G_dummy)

                # process cycles and correct order (according to 'preferential removal')
                if cyclecheck == True and ordercycles == True:
                    cycledict, degreedict = cycle_checks(G_dummy)

                # dijkstra
                ssource = nx.single_source_dijkstra(
                    G_dummy, source=root, weight="weight"
                )
                # result of dijkstra algorithm is sorted
                sortedssource = {
                    k: v
                    for k, v in sorted(ssource[0].items(), key=lambda item: item[1])
                }

                # get keys of sorted dict
                sortedssource_edges = sortedssource.keys()

                sortedssource_edges_list = list(sortedssource_edges)
                # sorted list contains the mutation route
                nodes = sortedssource_edges_list

                # order has to be reversed - the most distant atom is the first to be removed
                nodes.reverse()

                # sort nodes according to degree, cycle participation and removal order
                if cyclecheck == True and ordercycles == True:
                    nodes = change_route_cycles(
                        nodes, cycledict, degreedict, sortedssource, G
                    )

                logger.info("Final mutation route:")
                logger.info(nodes)
                ordered_LJ_mutations.append(nodes)

    return ordered_LJ_mutations


def _calculate_order_of_LJ_mutations_new_iter(
    connected_dummy_regions: list,
    match_terminal_atoms: dict,
    G: nx.Graph,
    cyclecheck=True,
    ordercheck=True,
) -> list:
    """
    bfs/djikstra-algorithm applied iteratively, i.e. after each removal of an atom
    """

    ordered_LJ_mutations = []

    for real_atom in match_terminal_atoms:
        for dummy_atom in match_terminal_atoms[real_atom]:
            for connected_dummy_region in connected_dummy_regions:
                # stop at connected dummy region with specific dummy_atom in it
                if dummy_atom not in connected_dummy_region:
                    continue

                G_dummy = G.copy()
                # delete all nodes not in dummy region
                remove_nodes = [
                    node for node in G.nodes() if node not in connected_dummy_region
                ]
                for remove_node in remove_nodes:
                    G_dummy.remove_node(remove_node)

                # root is the dummy atom that connects the real region with the dummy region
                root = dummy_atom

                final_order = []

                removeG = nx.Graph()
                removearray = []
                while len(G_dummy.nodes()) > 0:
                    # update weights according to already removed nodes
                    if ordercheck == True:
                        G_dummy = order_checks_nx(G_dummy, removearray, G)

                    G_origweights = G_dummy.copy()

                    # update weights according to cycle participation
                    if cyclecheck == True:
                        G_dummy = cycle_checks_nx(G_dummy)

                    # dijkstra
                    ssource = nx.single_source_dijkstra(
                        G_dummy, source=root, weight="weight"
                    )

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

                logger.info("Final mutation route:")
                logger.info(nodes)

                # no reverse (already correct order, starting with the node with greatest distance from root)
                ordered_LJ_mutations.append(nodes)

    return ordered_LJ_mutations


def _calculate_order_of_LJ_mutations_new_iter_change(
    connected_dummy_regions: list,
    match_terminal_atoms: dict,
    G: nx.Graph,
    cyclecheck=True,
    ordercheck=True,
) -> list:
    """
    works iteratively, i.e. after each removal of an atom, algorithm is chosen depending on current state
    """

    ordered_LJ_mutations = []

    initial_cdict, initial_degreedict = cycle_checks(G)

    cycles = nx.cycle_basis(G)

    for real_atom in match_terminal_atoms:
        for dummy_atom in match_terminal_atoms[real_atom]:
            for connected_dummy_region in connected_dummy_regions:
                # stop at connected dummy region with specific dummy_atom in it
                if dummy_atom not in connected_dummy_region:
                    continue

                G_dummy = G.copy()

                G_orig = G.copy()
                # delete all nodes not in dummy region
                remove_nodes = [
                    node for node in G.nodes() if node not in connected_dummy_region
                ]
                for remove_node in remove_nodes:
                    G_dummy.remove_node(remove_node)

                # root is the dummy atom that connects the real region with the dummy region
                root = dummy_atom

                # these variables will be changed accordingly to the current position (i.e. the last removed node) in the graph
                dfs_step = False
                cycle_step = False
                cycle_step_initialized = False
                illegal_root_cycle = False

                final_order = []

                removeG = nx.Graph()
                removearray = []

                while len(G_dummy.nodes()) > 0:

                    # update weights according to already removed nodes
                    if ordercheck == True:
                        G_dummy = order_checks_nx(G_dummy, removearray, G)

                    G_origweights = G_dummy.copy()

                    # update weights according to cycle participation
                    if cyclecheck == True:
                        G_dummy = cycle_checks_nx(G_dummy)

                    # if the last removed node is neither part of a longer chain nor of a cycle, an usual dijkstra step is carried out
                    if dfs_step == False and cycle_step == False:
                        # logger.info("dijkstra step")
                        # dijkstra
                        ssource = nx.single_source_dijkstra(
                            G_dummy, source=root, weight="weight"
                        )

                        sortedssource = {
                            k: v
                            for k, v in sorted(
                                ssource[0].items(), key=lambda item: item[1]
                            )
                        }
                        max_node = max(ssource[0], key=ssource[0].get)

                        current_degree = G_dummy.degree(max_node)

                        neighbors = [n for n in G_dummy.neighbors(max_node)]

                        # the current neighbor is determined
                        if current_degree > 0:
                            current_neighbor = neighbors[0]
                        else:
                            current_neighbor = None

                        # the type of the next iteration is determined
                        if current_neighbor != None:
                            if (
                                initial_cdict[max_node] == 0
                                and current_degree == 1
                                and G_dummy.degree(current_neighbor) == 2
                            ):
                                dfs_step = True

                            else:
                                illegal_root_cycle = False

                                for array in cycles:
                                    if max_node in array and root in array:
                                        illegal_root_cycle = True

                                if (
                                    initial_cdict[current_neighbor] > 0
                                    and illegal_root_cycle == False
                                ):

                                    cycle_step = True
                                    cycle_step_initialized = False

                                    cyclepath = sortedssource
                                    dfs_step = False

                            if current_neighbor == root:
                                dfs_step = False

                    # if the current position is within a cycle, a cycle_step is carried out
                    # neighbours within the cycle are determined and, if possible, removed
                    # it is not possible if the cycle atoms still have other bonds -> switch to dijkstra
                    elif cycle_step == True:
                        # logger.info("cycle step")
                        if cycle_step_initialized == False:

                            # get neighbors in cycle
                            for array in cycles:
                                if current_neighbor in array:
                                    currentcyclenodes = array

                            if max_node in currentcyclenodes:
                                currentcyclenodes.remove(max_node)

                            cyclepath2 = []
                            for el in cyclepath:
                                if el in currentcyclenodes:
                                    cyclepath2.append(el)

                            # check current state (cycle participation, degree)
                            current_cdict, current_degreedict = cycle_checks(G_dummy)

                            cyclepath_final = []

                            # check if the common core is correct, i.e. root does not participate in cycle
                            illegal_root_cycle = False
                            for node in currentcyclenodes:
                                if node not in connected_dummy_region:
                                    illegal_root_cycle = True
                                    dfs_step = False
                                    cycle_step = False
                                    cycle_step_initialized = False
                                    G_dummy = G_origweights.copy()
                                    continue

                            if len(cyclepath2) == 0:
                                dfs_step = False
                                cycle_step = False
                                cycle_step_initialized = False
                                G_dummy = G_origweights.copy()
                                continue

                            min_degree = current_degreedict[cyclepath2[0]]
                            for el in cyclepath2:
                                if current_degreedict[el] < min_degree:
                                    min_degree = current_degreedict[el]

                            # if degree is too high, switch to dijkstra
                            if min_degree > 1:
                                dfs_step = False
                                cycle_step = False
                                cycle_step_initialized = False
                                G_dummy = G_origweights.copy()
                                continue

                            for el in cyclepath2:
                                if current_degreedict[el] == min_degree:
                                    cyclepath_final.append(el)

                            cycle_step_initialized = True

                        problematic_neighbors = False
                        neighbors = [n for n in G_dummy.neighbors(cyclepath_final[-1])]

                        if len(neighbors) > 0:
                            for neighbor in neighbors:

                                if (
                                    initial_cdict[neighbor] == 0
                                    or initial_cdict[max_node] == 0
                                ):
                                    problematic_neighbors = True

                                    continue
                                # only for pathological cases (pdbbind l1a4k)
                                if (
                                    initial_cdict[max_node] > 2
                                    or initial_cdict[cyclepath_final[-1]] > 2
                                ) and G.degree(neighbor) < 3:
                                    problematic_neighbors = True

                                    continue

                        for elem in cyclepath:
                            if G.degree(elem) > 3:
                                problematic_neighbors = True

                        if problematic_neighbors == True:

                            dfs_step = False
                            cycle_step = False
                            cycle_step_initialized = False

                            G_dummy = G_origweights.copy()
                            continue

                        # only for pathological cases (root pertains to a cycle)
                        if root in cyclepath_final:
                            dfs_step = False
                            cycle_step = False
                            cycle_step_initialized = False
                            G_dummy = G_origweights.copy()
                            continue

                        illegal_root_cycle = False

                        for array in cycles:
                            if (
                                current_neighbor in array
                                or max_node in array
                                or cyclepath_final[-1] in array
                            ) and root in array:

                                illegal_root_cycle = True
                                dfs_step = False
                                cycle_step = False
                                cycle_step_initialized = False
                                continue
                        if illegal_root_cycle == True:

                            G_dummy = G_origweights.copy()
                            continue

                        max_node = cyclepath_final.pop()

                        if len(cyclepath_final) == 0:
                            dfs_step = False
                            cycle_step = False
                            cycle_step_initialized = False

                    # if neither the conditions for a cycle_step are met nor the dijkstra algorithm has to be applied -> the current position is within a chain
                    # a dfs step is carried out
                    else:
                        # logger.info("dfs step")
                        max_node = current_neighbor

                        current_degree = G_dummy.degree(max_node)

                        if current_degree > 0:
                            current_neighbor = [n for n in G_dummy.neighbors(max_node)][
                                0
                            ]
                        else:
                            current_neighbor = None

                        if current_neighbor != None:
                            if (
                                initial_cdict[max_node] > 0
                                or current_degree != 1
                                or G_dummy.degree(current_neighbor) != 2
                            ):
                                dfs_step = False
                            current_cdict, current_degreedict = cycle_checks(G_dummy)
                            if current_cdict[max_node] > 0:
                                ssource = nx.single_source_dijkstra(
                                    G_dummy, source=root, weight="weight"
                                )
                                sortedssource = {
                                    k: v
                                    for k, v in sorted(
                                        ssource[0].items(), key=lambda item: item[1]
                                    )
                                }
                                cyclepath = sortedssource
                                cycle_step = True
                        if current_neighbor == root:
                            dfs_step = False

                    if max_node == root and len(G_dummy.nodes()) > 1:
                        illegal_root_cycle = True
                        dfs_step = False
                        cycle_step = False
                        cycle_step_initialized = False
                        G_dummy = G_origweights.copy()

                        continue

                    # the determined node is added to the final array
                    final_order.extend([max_node])

                    G_dummy = G_origweights.copy()

                    # remove G_dummy
                    G_dummy.remove_node(max_node)

                    # add to removeG
                    removeG.add_node(max_node)
                    removearray.append(max_node)

                sortedssource_edges = final_order

                # sortedssource_edges_list already contains the nodes in right (reversed) order (including root)
                nodes = sortedssource_edges

                logger.info("Final mutation route:")
                logger.info(nodes)

                # no reverse (already correct order, starting with the node with greatest distance from root)
                ordered_LJ_mutations.append(nodes)

    return ordered_LJ_mutations


def cycle_checks_nx(G):
    """
    cycle processing, currently used in _calculate_order_of_LJ_mutations_new_iter
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
            G[el[0]][el[1]]["weight"] = G[el[0]][el[1]]["weight"] - cdict[i] * 5

    return G


def order_checks_nx(G, removearray, G_total):
    """
    preferential removal (atoms which neighbours already have been removed are removed earlier), currently used in _calculate_order_of_LJ_mutations_new_iter
    ------
    returns nx-graph-object with updated weights (according to removearray)
    """

    if len(removearray) > 0:
        lastremoved = removearray[len(removearray) - 1]

        edg = G_total.edges(lastremoved)

        edg_dummy = G.edges()

        for ed in edg:

            if ed[0] != lastremoved:
                connectednode = ed[0]
            else:
                connectednode = ed[1]

            # if node is connected to last removed node, its weight get a small increase
            if G.has_node(connectednode):

                connectededges = G.edges(connectednode)
                for el in connectededges:

                    G[el[0]][el[1]]["weight"] = G[el[0]][el[1]]["weight"] + 1

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
