from typing import Any, Dict, List, Tuple, Union

import pandas as pd
from networkx import (
    DiGraph,
    ancestors,
    edge_betweenness_centrality,
    is_directed_acyclic_graph,
    isolates,
    set_edge_attributes,
)
from networkx.algorithms import bipartite
from utils import graph_hash_utils
import logging
logger = logging.getLogger(__name__)


def generate_combination_graphs(G: DiGraph, method: str = "ebc", max_nr: int = 0) -> List[DiGraph]:
    """
    Generate combination graphs from a directed graph.

    Args:
        G (DiGraph): A directed graph.
        method (str): The method to use for generating combination graphs. Valid values: ["ebc", "all"].
        max_nr (int): The maximum number of combination graphs to generate.

    Returns:
        List[DiGraph]: A list of combination graphs.

    Raises:
        ValueError: If an invalid combination graphs generating method is provided.
    """
    logger.debug(f"Generating combination graphs using method: {method}")
    if method == "all":
        combination_graphs = generate_combination_graphs_all(G)
    elif method == "ebc":
        combination_graphs = select_edges_by_edge_betweenness_centrality(G, max_nr)
    else:
        # TODO: Replace with a custom exception
        raise ValueError("[ERROR] Invalid combination graphs generating method was provided. Valid values: [ebc|all].")

    return combination_graphs


def generate_combination_graphs_all(G: DiGraph) -> List[DiGraph]:
    """
    Generate all possible combination graphs from a directed graph.

    Args:
        G (DiGraph): A directed graph.

    Returns:
        List[DiGraph]: A list of combination graphs.
    """
    # 1. Assemble edges from those that should always be included (constant_edges) and from those that
    # belong to substances of >1 in_degree .
    # From the latter type of edges only one edge is kept, for this all combinations are enumerated, and these will
    # be merged with the constant_edges.
    #
    # 2. Then, an induced subgraph is created with the help of the merged edge list assembled above.
    #

    ###
    # 1
    ###

    # Identify edges that belong to substances that have out_degree = 1 property
    const_out_edges_substances = []
    const_out_edges_substances = get_all_out_edges_of_substances(G)

    # Identify edges that belong to substances that have out_degree > 1 property
    const_edges_substance_indegree_one = []
    const_edges_substance_indegree_one = get_edges_of_substances_indegree_one(G)

    # Identify edges that belong to substances that have out_degree > 1 property
    # These are the substances that are reactants/reagents of multiple reactions

    # edges_substance_outdegree_multiple = []
    # edges_substance_outdegree_multiple = get_edges_of_substances_outdegree_multiple (G)

    # Identify edges that belong to substances that have out_degree > 1 property
    # These are the substances that are products of multiple reactions

    edges_substance_indegree_multiple = []
    edges_substance_indegree_multiple = get_edges_of_substances_indegree_multiple(G)

    ###
    # 2
    ###

    constant_edges = []
    constant_edges = const_out_edges_substances
    constant_edges.extend(const_edges_substance_indegree_one)

    select_edges = []
    select_edges = edges_substance_indegree_multiple
    # select_edges = edges_substance_outdegree_multiple
    # edges_substance_outdegree_multiple.extend (edges_substance_indegree_multiple)

    select_combinations = get_select_edge_combinations(select_edges)

    final_combinations = []

    for combination in select_combinations:
        combination.extend(constant_edges)
        final_combinations.append(combination)
        # write_to_log (final_combinations)

    combination_graphs = []

    for combination in final_combinations:
        C = G.edge_subgraph(combination).copy()
        combination_graphs.append(C)

    logger.debug(f"Generated {len(combination_graphs)} combination graphs.")

    return combination_graphs


def select_edges_by_edge_betweenness_centrality(C: DiGraph, nr_variants: int = 0) -> List[DiGraph]:
    """
    Select edges by edge betweenness centrality.

    Args:
        C (DiGraph): A directed graph.
        nr_variants (int): The number of combination graphs to generate.

    Returns:
        List[DiGraph]: A list of combination graphs.
    """
    logger.info("Starting select_edges_by_edge_betweenness_centrality")
    G = C.copy()
    combination_graphs = []
    orig_nr_variants = nr_variants

    if nr_variants < 1:
        idx = -1
    else:
        idx = 0

    edge_choices = get_edges_of_substances_indegree_multiple(G)
    nr_multi_indegree_substances = len(edge_choices)

    edge_betweenness = edge_betweenness_centrality(G, k=None, normalized=True, weight=None, seed=None)

    while nr_multi_indegree_substances > 0 and idx < nr_variants:

        const_out_edges_substances = get_all_out_edges_of_substances(G)
        const_edges_substance_indegree_one = get_edges_of_substances_indegree_one(G)
        constant_edges = const_out_edges_substances + const_edges_substance_indegree_one

        start_nodes, end_nodes, ebcs = [], [], []
        for ec in edge_choices:
            for e in ec:
                start_nodes.append(e[0])
                end_nodes.append(e[1])
                ebcs.append(edge_betweenness[e])

        df = pd.DataFrame({"start_node": start_nodes, "end_node": end_nodes, "ebc": ebcs})
        df = df.sort_values(by=["end_node", "ebc", "start_node"], ascending=[True, False, True])
        df = df.groupby(["end_node"], as_index=False).agg("first")

        selected_edges = list(zip(df["start_node"], df["end_node"]))
        edges_to_delete = selected_edges.copy()
        selected_edges.extend(constant_edges)

        CG = G.edge_subgraph(selected_edges).copy()
        combination_graphs.append(CG)

        G.remove_edges_from(edges_to_delete)

        edge_choices = get_edges_of_substances_indegree_multiple(G)
        nr_multi_indegree_substances = len(edge_choices)

        edge_betweenness = edge_betweenness_centrality(G, k=None, normalized=True, weight=None, seed=None)

        if orig_nr_variants < 1:
            idx = -1
        else:
            idx += 1

    combination_graphs.append(G.copy())

    return combination_graphs


def get_select_edge_combinations(edges):
    # "This code is contributed by mohit kumar"

    # write_to_log (arr)

    # "number of arrays"
    n = len(edges)

    # "to keep track of next element
    # in each of the n arrays"
    indices = [0 for i in range(n)]

    results = []
    tmp = []  # noqa: F841

    while 1:
        combinaton = []

        # "write_to_log current combination"
        for i in range(n):
            combinaton.append(edges[i][indices[i]])

        # write_to_log()

        results.append(combinaton)
        # idx += 1

        # results.append(tmp)

        # "find the rightmost array that has more
        # elements left after the current element
        # in that array"
        next = n - 1
        while next >= 0 and (indices[next] + 1 >= len(edges[next])):
            next -= 1

        # "no such array is found so no more
        # combinations left"
        if next < 0:
            return results

        # "if found move to next element in that
        # array"
        indices[next] += 1

        # "for all arrays to the right of this
        # array current index again points to
        # first element"
        for i in range(next + 1, n):
            indices[i] = 0

    return None


def get_zero_yield_reactions(G: DiGraph) -> List:
    """
    Identify reaction node IDs with a predicted yield of 0.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List: A list of reaction node IDs where yield_predicted is 0.
    """
    return [n for n in G.nodes if G.nodes[n]["node_type"].lower() == "reaction" and G.nodes[n]["yield_predicted"] == 0]


def get_incompatible_reactions(G: DiGraph, ambiguous_rxns: List[str], multiproduct_rxns: List[str]) -> List:
    """
    Identify reactions node IDs that are incompatible with the target molecule.

    Args:
        G (DiGraph): A NetworkX directed graph.
        ambiguous_rxns (List[str]): A list of reaction IDs that are ambiguous.
        multiproduct_rxns (List[str]): A list of reaction IDs that produce multiple products.

    Returns:
        List: A list of reaction node IDs that are incompatible with the target molecule.
    """
    return ambiguous_rxns + multiproduct_rxns


def get_parent_nodes_of_starting_materials(G: DiGraph) -> List:
    """
    Identify the parent nodes of starting materials in the graph.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List: A list of node IDs representing the parent nodes of starting materials.
    """
    parent_nodes = {
        predecessor
        for node, attrs in G.nodes(data=True)
        if attrs.get("node_type").lower() == "substance" and attrs.get("srole") == "sm"
        for predecessor in G.pred[node]
    }
    return list(parent_nodes)


def get_childrens_of_target_molecule(G: DiGraph) -> List:
    """
    Identify the children of the target molecule in the graph.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List: A list of node IDs representing the children of the target molecule.
    """
    return [
        child
        for node, attrs in G.nodes(data=True)
        if attrs.get("node_type").lower() == "substance" and attrs.get("srole") == "tm"
        for child in G.succ[node]
        if len(G.succ[node]) > 0
    ]


def get_multiproduct_rxns(G):
    # any reaction node that has out_degree > 1 is a multiproduct reaction
    return [n for n in G.nodes if G.nodes[n]["node_type"] == "reaction" and G.out_degree(n) > 1]


def get_ambiguous_rxns(G):
    # if any outgoing substance from reaction is same as incoming substance, then its an ambiguous reaction
    ambiguous_rxns = []
    for n in G.nodes:
        if G.nodes[n]["node_type"] == "reaction":
            succs = set(G.successors(n))
            preds = set(G.predecessors(n))
            for succ in succs:
                if succ in preds:
                    ambiguous_rxns.append(n)
                    break
    return ambiguous_rxns


def remove_incompatible_reactions(
    G: DiGraph,
    target_molecule_node_id: Union[str, int],
) -> DiGraph:
    """
    Remove incompatible reactions from the graph.

    Args:
        G (DiGraph): A NetworkX directed graph.
        target_molecule_node_id (Union[str, int]): The node ID of the target molecule.

    Returns:
        DiGraph: A directed graph with incompatible reactions removed.
    """
    ambiguous_rxns = get_ambiguous_rxns(G)
    multiproduct_rxns = get_multiproduct_rxns(G)
    incompatible_reactions = set(get_incompatible_reactions(G, ambiguous_rxns, multiproduct_rxns))
    incompatible_reactions.update(get_childrens_of_target_molecule(G))
    incompatible_reactions.update(get_parent_nodes_of_starting_materials(G))
    incompatible_reactions.update(get_zero_yield_reactions(G))

    # Remove incompatible reactions if they exist

    if incompatible_reactions:
        logger.debug(f"Removing {len(incompatible_reactions)} incompatible reactions from the Graph.")
        for nn in incompatible_reactions:
            G.remove_node(nn)

        # Keep ancestors of the target molecule in the graph
        if G.nodes:
            ancestors_tm = ancestors(G, target_molecule_node_id)
            ancestors_tm.add(target_molecule_node_id)
            G = G.subgraph(ancestors_tm).copy()

    # Remove any isolated nodes, this would be any remaining substances that are not ancestors of the target molecule
    G.remove_nodes_from(list(isolates(G)))

    return G


# TODO: Create the following function and replace repeated code
# Target molecule component
# ancestors_tm = ancestors(G, target_molecule_node_id)
# ancestors_tm.add(target_molecule_node_id)
# G = G.subgraph(ancestors_tm).copy()
# Find ancestors of TM and create node induced subgraph
# Def 12: In a synthesis graph G, let the target molecule component (TMC) be defined as the subgraph ,
# derived as follows. Let the target molecule and its ancestors constitute a node set.
# We use this node set to generate a node-induced subgraph of G. The component obtained such
# manner is S, i.e.: the TMC. Note, that S is an LCC.


def get_edges_of_substances_indegree_multiple(
    G: DiGraph,
) -> List[List[Tuple[Union[int, str], Union[int, str]]]]:
    """
    Return a list of edges that belong to substances with an indegree greater than 1.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List[List[Tuple[Union[int, str], Union[int, str]]]]: A list of edges that belong to substances with an indegree greater than 1.
    """
    edges = [[(source_node, n) for source_node in G.pred[n]] for n in G.nodes if G.nodes[n]["node_type"].lower() == "substance" and G.in_degree[n] > 1]
    return edges


def get_all_out_edges_of_substances(
    G: DiGraph,
) -> List[Tuple[Union[int, str], Union[int, str]]]:
    """
    Return all outgoing edges of substances in the graph.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List[Tuple[Union[int, str], Union[int, str]]]: A list of all outgoing edges of substances.
    """
    edges = [(n, target_node) for n, nbrs in G.adjacency() if G.nodes[n]["node_type"].lower() == "substance" for target_node in nbrs]
    return edges


def count_substances_with_multiple_incoming_edges(G: DiGraph) -> int:
    """
    Count the number of 'substance' nodes in the graph with more than one incoming edge.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        int: The count of 'substance' nodes with an indegree greater than 1.
    """
    return sum(1 for node in G.nodes if G.nodes[node]["node_type"].lower() == "substance" and G.in_degree[node] > 1)


def count_substances_with_multiple_outgoing_edges(G: DiGraph) -> int:
    """
    Count the number of 'substance' nodes in the graph with more than one outgoing edge.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        int: The count of 'substance' nodes with an outdegree greater than 1.
    """
    return sum(1 for node in G.nodes if G.nodes[node]["node_type"].lower() == "substance" and G.out_degree[node] > 1)


def get_edges_of_substances_indegree_one(
    G: DiGraph,
) -> List[Tuple[Union[int, str], Union[int, str]]]:
    edges = [(source_node, n) for n, preds in G.pred.items() if G.nodes[n]["node_type"].lower() == "substance" and G.in_degree[n] == 1 for source_node in preds]
    return edges


def get_child_nodes_of_leaf_intermediates(G: DiGraph) -> List:
    """
    Identify the child nodes of leaf intermediate substances in the graph.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List: A list of node IDs representing the child nodes of leaf intermediate substances.
    """
    leaf_intermediates = get_leaf_intermediate_substances(G)

    # Using a set comprehension to avoid duplicates
    reactions = {c for li in leaf_intermediates for c in G.succ[li]}

    return list(reactions)


def get_leaf_intermediate_substances(G: DiGraph) -> List[Union[int, str]]:
    """
    Identify leaf intermediate substances in the graph.

    Args:
        G (DiGraph): A NetworkX directed graph.

    Returns:
        List[Union[int, str]]: A list of node IDs representing leaf intermediate substances.
    """
    leaf_intermediate_substances = []

    for n in G.nodes:
        if G.nodes[n].get("node_type").lower() == "substance" and G.nodes[n].get("srole") == "im" and G.in_degree[n] == 0:
            leaf_intermediate_substances.append(n)

    return leaf_intermediate_substances


def find_viable_route(C: DiGraph, target_module_node_id: Union[str, int]) -> DiGraph:
    """
    Find a viable synthesis route in a combination graph.

    Args:
        C (DiGraph): A combination graph.
        target_module_node_id (Union[str, int]): The node ID of the target molecule.

    Returns:
        DiGraph: A directed graph representing a viable synthesis route.
    """
    # 1. In an iterative cycle,
    #    Marking reactions having leaf-intermediates as parent node.
    #    Then, identifying the ancestors of the target_molecule, then generating the node-induced subgraph
    #    of combination graph with the ancestors and the target molecule.
    #    Repeat these steps until no reactions (javing leaf-node as parent) are marked.
    # 2. If the node-induced graph consists of the target molecule only, then return an empty graph.
    # 3. If the node-induced graph is not a DAG return
    route_status = "Viable Route Candidate"
    G = C.copy()

    marked_reaction_nodes = []
    marked = True

    ###
    # 1
    ###
    while marked:
        # Mark reactions having leaf-intermediates as parent node and removes them
        # Child nodes of leaf intermediates are reactions that have a im substance as a reactant or reagent
        marked_reaction_nodes = get_child_nodes_of_leaf_intermediates(G)
        marked_reaction_nodes = list(set(marked_reaction_nodes))
        if len(marked_reaction_nodes) > 0:
            marked = True
            for nn in marked_reaction_nodes:
                G.remove_node(nn)
        else:
            marked = False

            # keep ancestors of tm in
        node_nr = len(G.nodes)

        # Keep ancestors of the target molecule in the graph
        # Removes all nodes that are not ancestors of the target molecule
        if node_nr > 0:
            ancestors_tm = ancestors(G, target_module_node_id)
            ancestors_tm.add(target_module_node_id)
            G = G.subgraph(ancestors_tm).copy()

    if target_module_node_id not in G.nodes:
        # TODO: Replace with a custom exception
        raise Exception("[ERROR] Something went wrong in identifying viable synthesis routes. aicp.py, fn: find_viable_route ()")

    ###
    # 3
    ###

    if len(G.nodes) == 1:
        logger.debug("[***] Filtered graph only contains target molecule.")
        route_status = "Route Candidate Only Contains Target Molecule"
        return DiGraph(), route_status

    ###
    # 4
    ###

    if not is_directed_acyclic_graph(G):
        logger.debug("[***] Filtered graph is not a DAG.")
        route_status = "Route Candidate is not a DAG"
        return DiGraph(), route_status

    return G, route_status


def project_to_reaction_network(G: DiGraph) -> DiGraph:
    """
    Project a bipartite graph to a reaction network based on the node type.

    Args:
        G (DiGraph): A directed bipartite graph.

    Returns:
        DiGraph: The projected digraph focusing on reaction nodes.
    """
    # Determine the two sets of nodes in the bipartite graph
    nodes_a, nodes_b = bipartite.sets(G)

    # Infer the type of nodes in the first set
    first_node_type = G.nodes[next(iter(nodes_a))]["node_type"].lower()

    # Choose the set with "reaction" node type for projection; default to nodes_a otherwise
    nodes_for_projection = nodes_a if first_node_type == "reaction" else nodes_b

    # Project the graph based on the chosen node set
    projected_graph = bipartite.projected_graph(G, nodes_for_projection)

    # Tag edges in the projected graph as "projected"
    set_edge_attributes(projected_graph, "projected", "edge_type")

    return projected_graph


def extract_reaction_nodes_from_nxgraph(G: DiGraph):
    """
    Extracts reaction nodes from a directed graph.

    Args:
        G (DiGraph): A directed graph containing nodes with 'node_type' and 'rxid' attributes.

    Returns:
        List[str]: A list of reaction IDs.
    """
    return [G.nodes[node]["rxid"] for node in G.nodes if G.nodes[node]["node_type"].lower() == "reaction"]


def grab_yield_predicted(G: DiGraph) -> Dict[str, float]:
    """
    Extracts predicted yields for reactions from a directed graph.

    Args:
        G (DiGraph): A directed graph containing nodes with 'node_type', 'rxid', and 'yield_predicted' attributes.

    Returns:
        Dict[str, float]: A dictionary mapping reaction IDs to their predicted yields. If 'yield_predicted' is
                          the string 'False', it is treated as 0.0; otherwise, it is converted to a float.
    """
    inchikey_to_yield_predicted = {}
    for node in G.nodes:
        if G.nodes[node]["node_type"].lower() == "reaction":
            rxid = G.nodes[node]["node_label"]
            yield_predicted = G.nodes[node].get("yield_predicted")

            # Convert 'yield_predicted' to appropriate numerical value
            if yield_predicted == "False":
                inchikey_to_yield_predicted[rxid] = 0.0
            else:
                # Ensure yield_predicted is converted to float if not 'False'
                inchikey_to_yield_predicted[rxid] = float(yield_predicted) if yield_predicted is not None else 0.0

    return inchikey_to_yield_predicted


def annotate_projected_route_with_yield_info(route: DiGraph, yields: Dict[str, float]) -> None:
    """
    Annotates a projected route graph with yield information.

    Args:
        route (DiGraph): A projected route graph.
        yields (Dict[str, float]): A dictionary mapping rxid to yield percentages.

    Returns:
        None: The projected route graph is annotated in place.
    """
    yield_factor = 0.01

    # Annotate each node in the projected route with its yield information
    for node in route.nodes():
        logger.debug(route.nodes[node])
        rxid = route.nodes[node]["node_label"]
        # Convert yield percentage to a decimal and store it
        route.nodes[node]["yield"] = yields[rxid] * yield_factor


def identify_individual_synthesis_routes(
    G: DiGraph,
    target_molecule_node_id: Any,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Identify individual synthesis routes from a synthesis graph.

    Args:
        G (DiGraph): A directed graph representing a synthesis graph.
        target_molecule_node_id (Any): The node ID of the target molecule.

    Returns:
        Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]: A list of dictionaries containing individual synthesis routes.
    """
    #
    # Identifying all viable synthesis routes from a synthesis graph, and return them as
    # directed-acyclic-graphs (DAGs).
    #
    # Assumption: All reactions have a single products.
    #
    # 1. Remove reactions nodes; children of the target molecule, ambiguous reactions, multi-product reactions,
    #    and parents of starting materials.
    #    Generate a node-induced subgraph from the modified graph, using the target molecule and its ancestors.
    # 2. Generate combination graphs.
    # 3. Find viable synthesis route in each combination graph.
    # 4. Deduplicate the identified individual synhesis routes by sorting the node labels alphanumerically,
    #    then concatenating them by a comma in to a string, which uniqually identifies the
    #    synthesis route (which must be a DAG at this stage) at hand.
    #
    ###
    # 1
    ###
    # This now happens at an earlier stage in the pipeline
    # Leaving code for reference
    # G = remove_incompatible_reactions(
    #     G,
    #     target_molecule_node_id,
    #     ambiguous_rxns,
    #     multiproduct_rxns
    # )

    ###
    # 2
    ###

    combination_graphs = generate_combination_graphs(G)
    logger.debug(f"Identified {len(combination_graphs)} combination graphs.")

    ###
    # 3
    ###

    # Initialize a list to hold route candidates.
    route_candidates = []

    # Loop through each combination graph with an index starting from 1.
    for idx, combination_graph in enumerate(combination_graphs, start=1):
        # Initialize a dictionary to store the current route candidate.
        route_candidate: dict = {"route_index": idx}

        # If the current combination graph has nodes, attempt to find a viable route.
        # Otherwise, assign an empty DiGraph.
        if len(combination_graph.nodes) > 0:
            logger.debug(f"Identifying viable routes for combination graph {idx}.")
            route_candidate["route_candidate"], route_candidate["route_status"] = find_viable_route(combination_graph, target_molecule_node_id)
        else:
            logger.debug(f"Combination graph {idx} is empty.")
            route_candidate["route_candidate"] = DiGraph()
            route_candidate["route_status"] = "Empty Combination Graph"
        route_candidate["method"] = "AICP"

        # Add the current route candidate to the list.
        route_candidates.append(route_candidate)

    logger.debug(f"Identified {len(route_candidates)} route candidates.")

    ###
    # 4
    ###

    # Initialize a set to track unique tree hashes for efficient lookup.
    deduplicated_hashes = set()

    # Initialize a list to hold the final, deduplicated route candidates.
    final_route_candidates = []

    # Loop through each route candidate to generate tree hashes and deduplicate simultaneously.
    for route_candidate in route_candidates:
        graph = route_candidate["route_candidate"].copy()
        # Check if the graph is not empty.
        if graph.nodes:
            # Generate a tree hash for the graph.
            tree_hash = graph_hash_utils.hash_graph(graph)
            # Deduplicate based on the tree hash.
            if tree_hash not in deduplicated_hashes:
                # Mark this tree hash as encountered.
                deduplicated_hashes.add(tree_hash)
                # Store the tree and its hash in the route candidate.
                route_candidate["tree_hash"] = tree_hash
                # Add the route candidate to the final list.
                final_route_candidates.append(route_candidate)

    logger.debug(f"Identified {len(final_route_candidates)} final route candidates.")

    return final_route_candidates, route_candidates, combination_graphs
