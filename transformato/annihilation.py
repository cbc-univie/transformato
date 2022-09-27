import logging
import parmed as pm

import networkx as nx


logger = logging.getLogger(__name__)


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
