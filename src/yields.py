from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
from networkx import DiGraph


def aggregate_yields(G: DiGraph) -> float:
    (G, A, A_bin, node_map, node_reverse_map) = synth_route_to_adjacency_matrix(G)

    # Initialize row vector of aggregated yields
    # - identify the terminal reaction Tr (of which the product is the target molecule, TM)
    rt = identify_terminal_rxn(G)
    if rt is None:
        raise ValueError("Terminal reaction not found in the graph.")
    rt_idx = node_reverse_map[rt]

    leaf_nodes = identify_leaf_nodes(G)

    leaf_node_indices = [node_reverse_map[n] for n in leaf_nodes]

    C = accumulate_scale(G, rt_idx, node_map, A, is_reference=False)
    C_ref = accumulate_scale(G, rt_idx, node_map, A_bin, is_reference=True)

    # compute ratio between the ideal and computed => report this as aggregated yield, done!

    accumulated_scale_of_leaves = report_accumulated_scale_of_leaves(leaf_node_indices, C)
    accumulated_scale_of_leaves_ref = report_accumulated_scale_of_leaves(leaf_node_indices, C_ref)

    aggregated_yield = float(100.0 * float(accumulated_scale_of_leaves_ref / accumulated_scale_of_leaves))

    return aggregated_yield


def synth_route_to_adjacency_matrix(G: nx.DiGraph) -> Tuple[nx.DiGraph, np.ndarray, np.ndarray, Dict[int, str], Dict[str, int]]:
    # reverse edges
    G = G.reverse()

    # keeping track of nodes
    node_map: Dict[int, str] = {}
    node_reverse_map: Dict[str, int] = {}

    for i, n in enumerate(G.nodes()):
        node_map[i] = n
        node_reverse_map[n] = i

    # computing 1/yield (accounting for Div by Zero and for non-sensical negative yields)
    nodes = list(G.nodes())
    for n in nodes:
        y = G.nodes[n].get("yield", None)
        G.nodes[n]["one_over_yield"] = 1 / y

    # assign edge weights by assuming the end node's yield value
    for start_node, end_node in G.edges():
        G.edges[start_node, end_node]["one_over_yield"] = G.nodes[end_node]["one_over_yield"]

    # convert to weighted adjacency matrix
    A = nx.adjacency_matrix(G, nodelist=nodes, dtype=float, weight="one_over_yield").toarray()

    # convert to binary adjacency matrix
    A_bin = nx.adjacency_matrix(G, nodelist=nodes, dtype=float).toarray()

    return G, A, A_bin, node_map, node_reverse_map


def identify_terminal_rxn(G: nx.DiGraph) -> Optional[str]:
    """
    Identifies the terminal reaction in a reverse synthesis route graph.

    The terminal reaction (root node) is the node without incoming edges.

    Parameters:
        G (nx.DiGraph): A directed graph representing the reversed synthesis route.

    Returns:
        Optional[str]: The root node (terminal reaction) or None if not found.
    """
    root_node_idx: Optional[str] = None

    for n in G.nodes():
        if G.in_degree(n) == 0:
            if root_node_idx is not None:
                raise ValueError(f"[ERROR] There should be only one root node in the projected network, " f"but multiple were found. Reaction rxid: {n}")
            root_node_idx = n

    return root_node_idx


def identify_leaf_nodes(G: nx.DiGraph) -> List[str]:
    """
    Identifies the leaf nodes in a reverse synthesis route graph.

    Leaf nodes are defined as nodes with no outgoing edges.

    Parameters:
        G (nx.DiGraph): A directed graph representing the reversed synthesis route.

    Returns:
        List[str]: A list of leaf node identifiers.

    Raises:
        ValueError: If no leaf nodes are found in the projected network.
    """
    leaf_nodes: List[str] = [n for n in G.nodes() if G.out_degree(n) == 0]

    if not leaf_nodes:
        raise ValueError("[ERROR] There is no leaf node in the projected network!")

    return leaf_nodes


def accumulate_scale(G: nx.DiGraph, rt_idx: int, node_map: Dict[int, str], A: np.ndarray, is_reference: bool = False) -> np.ndarray:
    """
    Computes the accumulated scale of reactions in a reversed projected reaction network.

    Parameters:
        G (nx.DiGraph): Input reversed projected reaction network.
        rt_idx (int): Index of the terminal reaction (synthesizing the target molecule).
        node_map (Dict[int, str]): Mapping of node indices to node IDs.
        A (np.ndarray): Input adjacency matrix (either weighted or binary, but float type).
        is_reference (bool, optional): Whether to compute based on yield (False, default)
                                       or an ideal 100% yield scenario (True).

    Returns:
        np.ndarray: The accumulated scale of reactions.
    """
    num_nodes = len(G.nodes())
    row_vector = np.zeros(num_nodes, dtype=float)

    if is_reference:
        row_vector[rt_idx] = 1.0
    else:
        rt = node_map[rt_idx]
        row_vector[rt_idx] = G.nodes[rt]["one_over_yield"]

    # Initialize accumulated scale
    C = row_vector.copy()

    # Compute the longest path in the DAG as a boundary for iterations
    max_path_length = nx.dag_longest_path_length(G)

    for _ in range(max_path_length):
        row_vector = np.matmul(row_vector, A)
        C += row_vector

    return C


def report_accumulated_scale_of_leaves(leaf_node_indices: List[int], row_vector: np.ndarray) -> float:
    """
    Reports the accumulated scale of leaf nodes in the reaction network.

    Parameters:
        leaf_node_indices (List[int]): Indices of the leaf nodes.
        row_vector (np.ndarray): The row vector representing the accumulated scales.

    Returns:
        float: The total accumulated scale of the leaf nodes.
    """
    return float(sum(row_vector[idx] for idx in leaf_node_indices))
