import uuid
from typing import Any, Dict, List, Tuple

from networkx import DiGraph
from rdkit import Chem
from syngps.models.models import Provenance, RouteAssemblyType
from app.logging_config import logger
from app.models import (
    NoPathsFoundInAskcosResponse,
    NoResultFoundInAskcosResponse,
    TreeSearchResponse,
)

ASKCOS_NODE_TYPES = ["reaction", "chemical"]
AICP_REACTION_NODE_TYPE = "Reaction"
AICP_SUBSTANCE_NODE_TYPE = "Substance"
AICP_PRODUCT_OF_EDGE_TYPE = "PRODUCT_OF"
AICP_REACTANT_OF_EDGE_TYPE = "REACTANT_OF"


def askcos_tree2synth_graph(tree: TreeSearchResponse, USE_RETRO_RXN_RENDERING: bool = False) -> DiGraph:
    """
    Converts an ASKCOS graph to an AICP SynthGraph. This utilizes the 'graph' property of a TreeSearchResponse object.
    """
    if not tree.result:
        raise ValueError("No result found in provided ASKCOS tree search response.")

    # Get original graph from ASKCOS
    logger.debug("Processing ASKCOS Graph. Debug logging will provide entire ASKCOS response.")
    askcos_graph = tree.result.graph

    if not askcos_graph:
        raise ValueError("No graph found in provided ASKCOS tree search response.")

    # Retrieve
    askcos_nodes = askcos_graph["nodes"]
    askcos_edges = askcos_graph["links"]

    # Check for empty nodes and edges
    if not askcos_nodes or not askcos_edges:
        raise ValueError("Invalid graph found in provided ASKCOS tree search response.")

    # Process nodes and edges
    synth_graph = process_askcos_nodes_and_edges(askcos_nodes, askcos_edges, USE_RETRO_RXN_RENDERING)

    # Assign synthesis roles
    synth_graph = assign_synth_roles(synth_graph)

    return synth_graph


def askcos_tree2synth_paths(tree: TreeSearchResponse, USE_RETRO_RXN_RENDERING: bool = False) -> List[Dict[str, Any]]:
    """
    Converts an ASKCOS graph to a list of AICP SynthRoutes. This utilizes the 'paths' property of a TreeSearchResponse object.
    """
    if not tree.result:
        raise NoResultFoundInAskcosResponse("No result found in provided ASKCOS tree search response.")

    # Get original paths from ASKCOS
    logger.debug("Processing ASKCOS Paths. Debug logging will provide entire ASKCOS response.")
    askcos_paths = tree.result.paths

    if not askcos_paths:
        raise NoPathsFoundInAskcosResponse("No paths found in provided ASKCOS tree search response.")

    # Process each path
    synth_routes = []
    idx = 0
    for path in askcos_paths:
        # Process nodes and edges
        askcos_nodes = path["nodes"]
        askcos_edges = path["edges"]
        synth_graph = process_askcos_nodes_and_edges(askcos_nodes, askcos_edges, USE_RETRO_RXN_RENDERING)

        # Assign synthesis roles
        synth_graph = assign_synth_roles(synth_graph)

        # Create a SynthRoute object
        synth_route = {"idx": idx, "path": synth_graph}
        idx += 1

        # Append to list of SynthRoutes
        synth_routes.append(synth_route)

    return synth_routes


# Provided by Gergely in original ASCKOS2AICP codebase.
# Mostly by ChatGPT after several iterations:
def swap_reaction_sections(rxsmiles):
    # Separate rxsmiles and extension
    if "|" in rxsmiles:
        # Raise error if rxsmiles only has one pipe
        if rxsmiles.count("|") != 2:
            raise ValueError("Invalid RXSMILES format")
        reaction, extension = rxsmiles.split("|", 1)

        # Split RXSMILES into sections

        sections = reaction.strip().split(">")
        if len(sections) != 3:
            raise ValueError("Invalid RXSMILES format")
        reactants, reagents, products = sections
        reactant_components = reactants.split(".")
        reagent_components = reagents.split(".")
        product_components = products.split(".")
        # Parse the 'f:' extension
        ext_type, ext_data = extension.split(":")
        if ext_type != "f":
            raise ValueError("Unsupported extension type")
        component_groups = [list(map(int, group.split("."))) for group in ext_data.strip("|").split(",")]
        # Total components in each section
        total_reactants = len(reactant_components)
        total_reagents = len(reagent_components)
        total_products = len(product_components)
        # New boundaries after the swap
        new_total_reactants = total_products
        new_total_products = total_reactants
        # Update f: groups for the swap
        new_f = []
        for group in component_groups:
            # Adjust indices based on the section
            if all(i < total_reactants for i in group):  # Reactants -> Products
                new_f.append(".".join(str(i + total_reagents + total_products) for i in group))
            elif all(total_reactants <= i < total_reactants + total_reagents for i in group):  # Reagents
                # Reagents stay in the middle but need to shift to accommodate new boundaries
                new_f.append(".".join(str(i - total_reactants + new_total_reactants) for i in group))
            elif all(i >= total_reactants + total_reagents for i in group):  # Products -> Reactants
                new_f.append(".".join(str(i - total_reagents - total_reactants) for i in group))
            else:
                raise ValueError(f"Group {group} crosses section boundaries, violating rules.")
        # Rebuild the swapped RXSMILES
        swapped_rxsmiles = f"{products}>{reagents}>{reactants} |f:{','.join(new_f)}|"

    else:

        # Split RXSMILES into sections

        sections = rxsmiles.strip().split(">")

        if len(sections) != 3:
            raise ValueError("Invalid RXSMILES format")

        reactants, reagents, products = sections
        reactant_components = reactants.split(".")
        reagent_components = reagents.split(".")
        product_components = products.split(".")

        # Rebuild the swapped RXSMILES
        swapped_rxsmiles = f"{products}>{reagents}>{reactants}"

    return swapped_rxsmiles


def convert_askcos_edge_to_synth_edge(
    source: str, target: str, node_type_map: Dict[str, str], edge_type: str, USE_RETRO_RXN_RENDERING: bool = False
) -> Dict[str, Any]:
    """
    Converts an edge from an ASKCOS graph to an AICP SynthGraph edge.

    :param source: Source node ID
    :param target: Target node ID
    :param node_type_map: Map of node IDs to node types
    :param edge_type: Type of the edge (e.g., "product_of" or "reactant_of")
    :param USE_RETRO_RXN_RENDERING: Whether to use retro rendering for reactions
    :return: Edge metadata dictionary
    """
    start_node = source
    end_node = target

    # Handle reversed edges for reaction-substance relationships (already explained in the function)
    if node_type_map[source] == AICP_SUBSTANCE_NODE_TYPE and node_type_map[target] == AICP_REACTION_NODE_TYPE:
        edge_type = AICP_PRODUCT_OF_EDGE_TYPE
        start_node = target
        end_node = source
    elif node_type_map[source] == AICP_REACTION_NODE_TYPE and node_type_map[target] == AICP_SUBSTANCE_NODE_TYPE:
        edge_type = AICP_REACTANT_OF_EDGE_TYPE
        start_node = target
        end_node = source
    else:
        raise ValueError(f"Invalid edge type found between nodes {source} and {target}")

    # Directly create the edge metadata dictionary without using objects
    edge_metadata = {
        "uuid": f"{edge_type}_{uuid.uuid4().hex}",
        "start_node": start_node,
        "end_node": end_node,
        "edge_label": f"{start_node}|{end_node}",
        "edge_type": edge_type,
        "is_predicted": True,
        "inchikey": "",
        "rxid": "",
        "provenance": Provenance(is_in_askcos=True),
        "route_assembly_type": RouteAssemblyType(is_predicted=True, is_evidence=False),
    }

    return edge_metadata


def process_askcos_nodes_and_edges(askcos_nodes: List[Dict[Any, Any]], askcos_edges: List[Dict[str, str]], USE_RETRO_RXN_RENDERING: bool = False) -> DiGraph:
    graph = DiGraph()
    node_type_map = {}

    # Process all nodes
    for node in askcos_nodes:
        node_metadata = convert_askcos_node_to_synth_node(node, USE_RETRO_RXN_RENDERING)
        node_type_map[node_metadata["node_id"]] = node_metadata["node_type"]
        graph.add_node(node_metadata["node_id"], **node_metadata)

    # Process all edges
    for edge in askcos_edges:
        source = edge["source"] if "source" in edge else edge["from"]
        target = edge["target"] if "target" in edge else edge["to"]

        # Ensure source and target nodes exist in the graph
        if source not in graph.nodes:
            raise ValueError(f"Source node {source} not found in graph.")
        if target not in graph.nodes:
            raise ValueError(f"Target node {target} not found in graph.")

        # Convert edge using the Edge class and edge metadata
        edge_metadata = convert_askcos_edge_to_synth_edge(source, target, node_type_map, edge_type="reactant_of")

        # Add edge to graph
        graph.add_edge(edge_metadata["start_node"], edge_metadata["end_node"], **edge_metadata)

    # Return final generated graph
    return graph


def convert_askcos_node_to_synth_node(askcos_node: Dict[Any, Any], USE_RETRO_RXN_RENDERING: bool = False) -> Dict[Any, Any]:
    """
    Converts a node from an ASKCOS tree to a node recognized in an AICP SynthGraph.

    :param askcos_node: A node from an ASKCOS tree.
    :param USE_RETRO_RXN_RENDERING: Flag to enable retro reaction rendering for reactions.
    :return: Dictionary representation of the converted node.
    """
    node_type = askcos_node["type"]
    node_id = askcos_node["id"]  # Added line to include node_id

    if node_type not in ASKCOS_NODE_TYPES:
        raise ValueError(f"Invalid node type found from ASKCOS: {node_type}. ASKCOS Node ID: {node_id}")

    # Handle reaction nodes
    if node_type == "reaction":
        rxid = askcos_node["id"]
        yield_score = askcos_node.get("yield_score", 0.0)
        yield_predicted = askcos_node.get("yield_predicted", 0.0)
        rxsmiles = askcos_node.get("smiles", rxid)

        if USE_RETRO_RXN_RENDERING:
            rxsmiles = swap_reaction_sections(rxsmiles)

        # Directly build the reaction dictionary
        reaction_dict = {
            "node_id": node_id,
            "node_label": rxid,
            "uuid": f"reaction_{uuid.uuid4().hex}",
            "yield_predicted": yield_predicted,
            "yield_score": yield_score,
            "is_predicted": True,
            "is_balanced": askcos_node.get("is_balanced", False),
            "rxid": rxid,
            "rxsmiles": rxsmiles,
            "node_type": AICP_REACTION_NODE_TYPE,
            "provenance": Provenance(is_in_askcos=True),
            "route_assembly_type": RouteAssemblyType(is_predicted=True, is_evidence=False),
        }

        return reaction_dict

    # Handle substance nodes
    smiles = askcos_node.get("smiles", node_id)
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchikey = Chem.MolToInchiKey(mol)
    except Exception:
        raise ValueError(f"RDKit problem with parsing SMILES {smiles} and/or generating InChI-Key. ASKCOS Node ID: {node_id}")

    # Directly build the substance dictionary
    substance_dict = {
        "node_id": node_id,
        "node_label": node_id,
        "uuid": f"substance_{uuid.uuid4().hex}",
        "inchikey": inchikey,
        "canonical_smiles": Chem.CanonSmiles(smiles),
        "srole": askcos_node.get("srole", ""),
        "is_predicted": True,
        "node_type": AICP_SUBSTANCE_NODE_TYPE,
        "provenance": Provenance(is_in_askcos=True),
        "route_assembly_type": RouteAssemblyType(is_predicted=True, is_evidence=False),
    }

    return substance_dict


def identify_target_molecule(graph: DiGraph) -> str:
    """
    Returns the node ID of the target molecule. Target molecule will be the node with 0 outgoing edges of node type 'Substance'.

    :param graph: A DiGraph object representing a synthesis graph.

    :return: The node ID of the target molecule.
    """
    target_molecule = None
    for node in graph.nodes(data=True):
        if node[1]["node_type"] == AICP_SUBSTANCE_NODE_TYPE and graph.out_degree(node[0]) == 0:
            target_molecule = node[0]
            break

    if target_molecule is None:
        raise ValueError("No target molecule found in the provided synthesis graph.")

    return target_molecule


def assign_synth_roles(graph: DiGraph) -> DiGraph:
    """
    Assigns synthesis roles to nodes in a synthesis graph generated from an ASKCOS tree.

    :param graph: A DiGraph object representing a synthesis graph.

    :return: A DiGraph object with synthesis roles assigned to nodes.
    """
    target_node = identify_target_molecule(graph)

    # Assign target molecule role
    graph.nodes[target_node]["srole"] = "tm"

    # Assign im and sm roles
    for node in graph.nodes:
        node_data = graph.nodes[node]  # Access node attributes

        # Skip any reaction nodes
        if node_data["node_type"] == AICP_REACTION_NODE_TYPE:
            continue

        # Skip target molecule
        if node_data["node_id"] == target_node:
            continue

        # If the node has 0 incoming edges, mark it as a sm (starting material), otherwise as an im (intermediate material)
        incoming_edges = graph.in_degree(node)
        if incoming_edges == 0:
            graph.nodes[node]["srole"] = "sm"
        else:
            graph.nodes[node]["srole"] = "im"

    return graph


def askcos_tree2synth_paths_with_graph(tree: TreeSearchResponse, USE_RETRO_RXN_RENDERING: bool = False) -> Tuple[DiGraph, List[Dict[str, Any]]]:
    """
    Converts an ASKCOS graph to a synth graph and a list of AICP paths.
    """
    if not tree.result:
        raise NoResultFoundInAskcosResponse("No result found in provided ASKCOS tree search response.")

    # Get original paths from ASKCOS
    askcos_paths = tree.result.paths

    if not askcos_paths:
        raise NoPathsFoundInAskcosResponse("No paths found in provided ASKCOS tree search response.")

    graph_nodes = []
    graph_edges = []
    graph_paths = []

    # Step 0: Process each path
    path_index = 0
    for path in askcos_paths:
        path_nodes = []
        path_edges = []
        for node in path["nodes"]:
            graph_nodes.append(node)
            path_nodes.append(node)
        for edge in path["edges"]:
            graph_edges.append(edge)
            path_edges.append(edge)
        graph_paths.append({"nodes": path_nodes, "edges": path_edges, "path_index": path_index})
        path_index += 1

    # Step 1: Merge nodes with same SMILES
    smiles_to_node = {}
    old_id_to_new_id = {}
    for node in graph_nodes:
        smiles = node["smiles"]
        if smiles not in smiles_to_node:
            smiles_to_node[smiles] = node
        old_id_to_new_id[node["id"]] = smiles_to_node[smiles]["id"]

    merged_nodes = list(smiles_to_node.values())

    # Step 2: Update edges with new node IDs, swapping nodes
    for edge in graph_edges:
        edge["from"] = old_id_to_new_id[edge["from"]]
        edge["to"] = old_id_to_new_id[edge["to"]]

    # Step 3: Update paths to use new IDs
    for path in graph_paths:
        # Map nodes to just the IDs
        path["nodes"] = [old_id_to_new_id[node["id"]] for node in path["nodes"]]

        # Map edges to just the labels
        path["edges"] = [edge["to"] + "|" + edge["from"] for edge in path["edges"]]

    # Step 4: Synthesize graph
    synth_graph = process_askcos_nodes_and_edges(merged_nodes, graph_edges, USE_RETRO_RXN_RENDERING)

    # Step 5: Assign synthesis roles
    synth_graph = assign_synth_roles(synth_graph)

    # Return
    return (synth_graph, graph_paths)
