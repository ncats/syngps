from typing import Dict, List, Tuple, Union

from networkx import DiGraph
from networkx.readwrite import json_graph
import numpy as np
import datetime
import json
import logging

logger = logging.getLogger(__name__)


def json_routes_to_graph_subgraphs(json_graph_obj: Dict, routes_list: List) -> Tuple[DiGraph, List]:
    """
    Parse a graph_obj representation of a graph, process nodes, edges, and routes,
    then return the subgraphs for each route.

    Args:
        graph_obj (dict): A graph_obj representing the graph.

    Returns:
        dict: A dictionary with route names as keys and corresponding subgraphs as values.
    """
    # Parse graph
    G = json_to_graph(json_graph_obj)

    # Parse routes and create subgraphs from routes
    subgraphs = create_subgraphs_from_routes(routes_list, G)

    return G, subgraphs


def create_subgraphs_from_routes(routes: List, G: DiGraph) -> List:
    """
    Create subgraphs based on the nodes defined in the routes, including the edges
    that connect them automatically.

    Args:
        routes (dict): A dictionary containing routes information from the JSON.
        G (DiGraph): The original graph.

    Returns:
        dict: A dictionary with route indices as keys and corresponding subgraphs as values.
    """
    from syngps.models.models import (  # avoid circular imports
        SynthRoute,
    )

    subgraphs = []

    # Process each subgraph route in the 'subgraphs' list
    for route in routes:
        # Treat route as RouteLabelsInput
        route_data = route

        # Extract node labels for this route
        route_node_labels = route_data.get("route_node_labels", [])

        # Create a subgraph for these nodes
        subgraph = G.subgraph(route_node_labels).copy()

        # Create a SynthRoute object
        synth_route = SynthRoute(
            route_candidate=subgraph,
            aggregated_yield=route_data.get("aggregated_yield"),
            predicted=route_data.get("predicted", False),
            source=route_data.get("source"),
            route_index=route_data.get("route_index"),
            route_status=route_data.get("route_status"),
            method=route_data.get("method"),
        )

        # Append the subgraph to the list
        subgraphs.append(synth_route)

    return subgraphs


def graph_to_json(G: DiGraph) -> Dict:
    """
    Convert a NetworkX DiGraph to a JSON representation.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.

    Returns:
        dict: A JSON dictionary representing the graph.
    """
    graph_json = {
        "nodes": [{"node_label": node, **G.nodes[node]} for node in G.nodes()],
        "edges": [{"start_node": u, "end_node": v, "edge_label": f"{u}_{v}", **G.edges[u, v]} for u, v in G.edges()],
    }
    # logger.debug(f"Graph JSON: {graph_json}")
    return graph_json


def graph_to_json_nodes(G: DiGraph) -> List:
    """
    Convert a NetworkX DiGraph to a JSON representation.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.

    Returns:
        dict: A JSON dictionary representing the graph.
    """
    return [_parse_node({"node_label": node, **G.nodes[node]}) for node in G.nodes()]


def graph_to_json_edges(G: DiGraph) -> List:
    """
    Convert a NetworkX DiGraph to a JSON representation.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.

    Returns:
        dict: A JSON dictionary representing the graph.
    """
    return [_parse_edge({"start_node": u, "end_node": v, "node_label": f"{u}_{v}", **G.edges[u, v]}) for u, v in G.edges()]


def _parse_node(node_data: dict) -> dict:
    """
    Parses raw node data into a `ReactionNode` or `SubstanceNode`,
    replacing specific keys with grouped attributes like `provenance`,
    `validation`, and `route_assembly_type`. Ensures `node_type` is lowercase
    and removes blacklisted keys.

    Returns the dictionary representation of the node.
    """
    from syngps.models.models import RouteAssemblyType  # Import RouteAssemblyType
    from syngps.models.models import (  # avoid circular imports
        Provenance,
        ReactionNode,
        SubstanceNode,
        Validation,
        YieldInfo,
    )

    # Lists for dynamic assignments
    predicted_keys = ["is_in_savi_130k"]
    evidence_keys = ["is_in_uspto_full", "is_in_aicp"]
    black_key_list = ["node_id"]  # Keys to be removed

    # Determine `is_predicted` and `is_evidence` based on the presence of specific keys
    is_predicted = any(node_data.get(key) for key in predicted_keys)
    is_evidence = any(node_data.get(key) for key in evidence_keys)

    # Remove the original keys that are no longer needed
    original_data = node_data.copy()
    for key in predicted_keys + evidence_keys + black_key_list:
        node_data.pop(key, None)

    # Ensure `node_type` is lowercase if it exists
    if "node_type" in node_data:
        node_data["node_type"] = node_data["node_type"].lower()

    # Group `is_predicted` and `is_evidence` under `route_assembly_type` using RouteAssemblyType model
    route_assembly_type = RouteAssemblyType(
        is_predicted=is_predicted, is_evidence=is_evidence)

    node: Union[SubstanceNode, ReactionNode, None] = None

    # Create the appropriate node class based on `node_type`
    if node_data["node_type"] == "substance":
        node_data["route_assembly_type"] = route_assembly_type
        node_data["provenance"] = Provenance(
            patents=original_data["patents"] if "patents" in original_data else [
            ],
            patent_paragraph_nums=original_data["patent_paragraph_nums"] if "patent_paragraph_nums" in original_data else [
            ],
            is_in_aicp=original_data["is_in_aicp"] if "is_in_aicp" in original_data else False,
            is_in_savi_130k=original_data["is_in_savi_130k"] if "is_in_savi_130k" in original_data else False,
            is_in_uspto_full=original_data["is_in_uspto_full"] if "is_in_uspto_full" in original_data else False,
        )
        node = SubstanceNode(**node_data)

    elif node_data["node_type"] == "reaction":
        # Extract correct fields from node_data and map them to respective models
        node_data["yield_info"] = YieldInfo(
            yield_predicted=round(node_data["yield_predicted"], 2),
            yield_score=node_data["yield_score"],
        )
        node_data["provenance"] = Provenance(
            patents=original_data["patents"] if "patents" in original_data else [
            ],
            patent_paragraph_nums=original_data["patent_paragraph_nums"] if "patent_paragraph_nums" in original_data else [
            ],
            is_in_aicp=original_data["is_in_aicp"] if "is_in_aicp" in original_data else False,
            is_in_savi_130k=original_data["is_in_savi_130k"] if "is_in_savi_130k" in original_data else False,
            is_in_uspto_full=original_data["is_in_uspto_full"] if "is_in_uspto_full" in original_data else False,
        )
        if "is_rxname_recognized" not in node_data:
            node_data["is_rxname_recognized"] = False
        node_data["validation"] = Validation(
            is_balanced=node_data["is_balanced"],
            is_rxname_recognized=node_data["is_rxname_recognized"],
        )
        # Assign route_assembly_type
        node_data["route_assembly_type"] = route_assembly_type

        node = ReactionNode(**node_data)

    else:
        raise ValueError(f"Unknown node type: {node_data['node_type']}")

    # Return the dictionary representation of the instantiated object
    return node.dict()


def _parse_edge(edge_data: dict) -> dict:
    """
    Parses raw edge data, grouping attributes into `provenance`, `validation`,
    and `route_assembly_type`, and returns the dictionary representation of an Edge object.
    """
    from syngps.models.models import (  # avoid circular imports
        Edge,
        Provenance,
        RouteAssemblyType,
    )

    predicted_keys = ["is_in_savi_130k"]
    evidence_keys = ["is_in_uspto_full", "is_in_aicp"]

    # Determine `is_predicted` and `is_evidence`
    is_predicted = any(edge_data.get(key)
                       for key in predicted_keys) or edge_data.get("is_predicted", False)
    is_evidence = any(edge_data.get(key) for key in evidence_keys)

    # Ensure that `provenance` exists; default to `is_in_aicp=False`
    if "provenance" not in edge_data:
        edge_data["provenance"] = Provenance(
            is_in_aicp=edge_data.get("is_in_aicp", False),
            is_in_uspto_full=edge_data.get("is_in_uspto_full", False),
            is_in_savi_130k=edge_data.get("is_in_savi_130k", False),
        )

    # Blacklisted keys to remove
    black_key_list = ["start_node_id", "end_node_id", "node_label"]

    # Remove unnecessary keys
    for key in black_key_list + predicted_keys + evidence_keys:
        edge_data.pop(key, None)

    # Assign status information
    edge_data["route_assembly_type"] = RouteAssemblyType(
        is_predicted=is_predicted, is_evidence=is_evidence)

    # Create the Edge object from the processed data
    edge = Edge(**edge_data)

    # Return the dictionary representation of the instantiated Edge object
    return edge.dict()


def json_to_graph(graph_json: Dict) -> DiGraph:
    """
    Parse a JSON representation of a graph and return a NetworkX DiGraph.

    Args:
        graph_json (dict): A JSON dictionary representing the graph.

    Returns:
        DiGraph: A NetworkX DiGraph instance.
    """
    G = DiGraph()

    # Process nodes and edges
    process_nodes(G, graph_json["nodes"])
    process_edges(G, graph_json["edges"])

    return G


def add_node(G: DiGraph, node_id, node_metadata: dict) -> None:
    """
    Add a node with the specified ID and metadata to a directed graph.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.
        node_id (str): The ID of the node to be added.
        node_metadata (dict): A dictionary of node attributes.
    """
    node_metadata["node_id"] = node_id  # Ensure node_id is included in metadata
    G.add_node(node_id, **node_metadata)


def validate_node(node: Dict, valid_node_types: List) -> bool:
    """
    Validate if the node type is in the list of valid node types.

    Args:
        node (Dict): A dictionary representing a node.
        valid_node_types (List): A list of valid node types.

    Returns:
        bool: True if the node is valid, False otherwise.
    """
    if hasattr(node, "model_dump"):
        n = node.model_dump()
    elif isinstance(node, dict):
        n = node
    else:
        raise TypeError("Edge must be a Pydantic model or a dictionary")
    node_type = n.get("node_type")
    if node_type is None:
        raise ValueError("Node type is missing")
    node_type = node_type.lower()
    if node_type not in valid_node_types:
        raise ValueError(f"Invalid node type: {node_type}")
    return True


def process_nodes(G: DiGraph, nodes: List) -> DiGraph:
    """
    Process a list of nodes and add them to a directed graph.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.
        nodes (List): A list of dictionaries, each representing a node.

    Returns:
        DiGraph: The modified DiGraph with added nodes.
    """
    valid_node_types = ["reaction", "substance"]

    for node in nodes:
        if validate_node(node, valid_node_types):
            # Get node metadata
            if hasattr(node, "model_dump"):
                n = node.model_dump()
            elif isinstance(node, dict):
                n = node
            else:
                raise ValueError("Invalid node structure type")

            node_metadata = {"uuid": n.get("uuid"), "node_type": n.get(
                "node_type").lower(), "node_label": n.get("node_label")}

            if n.get("node_type").lower() == "reaction":
                node_metadata.update(
                    {
                        "rxid": n.get("rxid", None),
                        "rxsmiles": n.get("rxsmiles", None),
                        "name": n.get("node_label"),
                        "yield_score": n.get("yield_score", 0.0),
                        "yield_predicted": round(n.get("yield_predicted", 0.0), 2),
                        "yield_info": n.get("yield_info", {}),
                        "is_balanced": n.get("is_balanced", False),
                        "is_in_aicp": n.get("is_in_aicp", False),
                        "is_in_savi_130k": n.get("is_in_savi_130k", False),
                        "is_in_uspto_full": n.get("is_in_uspto_full", False),
                    }
                )
                add_node(G, n.get("node_label"), node_metadata)
            else:  # Substance node
                node_metadata.update(
                    {
                        "inchikey": n.get("inchikey", None),
                        "canonical_smiles": n.get("canonical_smiles", None),
                        "name": n.get("node_label"),
                        "is_in_aicp": n.get("is_in_aicp", False),
                        "is_in_savi_130k": n.get("is_in_savi_130k", False),
                        "is_in_uspto_full": n.get("is_in_uspto_full", False),
                        "srole": n.get("srole", None),
                    }
                )
                add_node(G, n.get("node_label"), node_metadata)

    return G


def add_edge(G: DiGraph, start_node, end_node, edge_metadata: dict) -> None:
    """
    Add an edge with metadata between two nodes in a directed graph.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.
        start_node (str): The starting node of the edge.
        end_node (str): The ending node of the edge.
        edge_metadata (dict): A dictionary of edge attributes.
    """
    G.add_edge(start_node, end_node, **edge_metadata)


def validate_edge(edge: dict, valid_edge_types: list) -> bool:
    """
    Validate if the node type is in the list of valid node types.

    Args:
        edge (dict): A dictionary representing an edge.
        valid_edge_types (list): A list of valid edge types.

    Returns:
        bool: True if the edge is valid, False otherwise.
    """
    if hasattr(edge, "model_dump"):
        e = edge.model_dump()
    elif isinstance(edge, dict):
        e = edge
    else:
        raise TypeError("Edge must be a Pydantic model or a dictionary")
    edge_type = e.get("edge_type").lower()
    if edge_type not in valid_edge_types:
        raise ValueError(f"Invalid edge type: {edge_type}")
    return True


def process_edges(G: DiGraph, edges: List) -> None:
    """
    Process a list of edges and add them to a directed graph.

    Args:
        G (DiGraph): A NetworkX DiGraph instance.
        edges (list): A list of dictionaries, each representing an edge.

    Returns:
        DiGraph: The modified DiGraph with added edges.
    """
    valid_edge_types = ["reactant_of", "reagent_of", "product_of"]
    for edge in edges:
        if validate_edge(edge, valid_edge_types):
            if hasattr(edge, "model_dump"):
                e = edge.model_dump()
            elif isinstance(edge, dict):
                e = edge
            else:
                raise TypeError(
                    "Edge must be a Pydantic model or a dictionary")
            edge_metadata = {}
            start_node = e["start_node"]
            end_node = e["end_node"]
            edge_metadata["edge_type"] = e["edge_type"]
            edge_metadata["uuid"] = e["uuid"]
            add_edge(G, start_node, end_node, edge_metadata)
    return G


def _make_json_safe(value):
    """
    Recursively converts a value to a JSON-serializable form.
    """
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    elif isinstance(value, (datetime.datetime, datetime.date)):
        return value.isoformat()
    elif isinstance(value, (np.integer, np.floating)):
        return value.item()
    elif isinstance(value, np.ndarray):
        return value.tolist()
    elif isinstance(value, dict):
        return {str(k): _make_json_safe(v) for k, v in value.items()}
    elif isinstance(value, (list, tuple, set)):
        return [_make_json_safe(v) for v in value]
    else:
        # Fallback for unknown types
        return str(value)


def safe_node_link_data(graph: DiGraph) -> dict:
    """
    Converts a NetworkX DiGraph to a JSON-serializable dict using node-link format.
    Handles datetime, NumPy types, and other common serialization issues.
    """
    if not isinstance(graph, DiGraph):
        raise TypeError("Graph must be a networkx.DiGraph")

    data = json_graph.node_link_data(graph)

    def safe_convert_attributes(obj):
        if isinstance(obj, dict):
            return {k: _make_json_safe(v) for k, v in obj.items()}
        return obj

    # Sanitize all node and edge attributes
    data["nodes"] = [safe_convert_attributes(
        node) for node in data.get("nodes", [])]
    data["links"] = [safe_convert_attributes(
        link) for link in data.get("links", [])]
    data["graph"] = safe_convert_attributes(data.get("graph", {}))

    return data
