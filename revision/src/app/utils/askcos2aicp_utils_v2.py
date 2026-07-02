import uuid
from typing import Any, Dict, List, Tuple
from app.utils.askcos2aicp_utils import (
    assign_synth_roles,
    swap_reaction_sections,
    convert_askcos_edge_to_synth_edge,
)
from networkx import DiGraph
from rdkit import Chem
from syngps.models.models import Provenance, RouteAssemblyType
from app.models import (
    NoPathsFoundInAskcosResponse,
)

# UDS (Universal Data Structure) parsing functions for the new ASKCOS response format.
# Reference: https://askcos-docs.mit.edu/guide/11-References/11.2-Universal-Data-Structure-(UDS).html
#
# UDS structure:
#   uds["node_dict"]          - {smiles: node_metadata} for all nodes in the search graph
#   uds["graph"]              - [{source: smiles, target: smiles}, ...] (product → reaction edges)
#   uds["uuid2smiles"]        - {uuid: smiles} mapping used by pathway edges
#   uds["pathways"]           - [[{source: uuid, target: uuid}, ...], ...] per-pathway edge lists
#   uds["pathways_properties"]- [{depth, precursor_cost, score, ...}, ...] per pathway
#
# Node IDs in UDS are SMILES strings (for graph/node_dict) or UUIDs (for pathways).
# Edges in pathways cover both directions:
#   - chemical_product → reaction  (substance is product of that reaction)
#   - reaction → reactant_chemical (reaction consumes that substance)

ASKCOS_NODE_TYPES = ["reaction", "chemical"]
AICP_REACTION_NODE_TYPE = "Reaction"
AICP_SUBSTANCE_NODE_TYPE = "Substance"
AICP_PRODUCT_OF_EDGE_TYPE = "PRODUCT_OF"
AICP_REACTANT_OF_EDGE_TYPE = "REACTANT_OF"

def uds_node_to_synth_node(smiles: str, node_data: Dict[Any, Any], USE_RETRO_RXN_RENDERING: bool = False) -> Dict[Any, Any]:
    """
    Converts a UDS node_dict entry to an AICP SynthGraph node dictionary.

    :param smiles: SMILES string used as the node ID in UDS.
    :param node_data: Metadata dict for the node from uds["node_dict"].
    :param USE_RETRO_RXN_RENDERING: If True, swaps reaction sections for retro rendering.
    :return: Dictionary representation of the converted node.
    """
    node_type = node_data.get("type")

    if node_type not in ASKCOS_NODE_TYPES:
        raise ValueError(f"Invalid node type '{node_type}' for SMILES: {smiles}")

    if node_type == "reaction":
        rxsmiles = node_data.get("smiles", smiles)
        if USE_RETRO_RXN_RENDERING:
            rxsmiles = swap_reaction_sections(rxsmiles)

        return {
            "node_id": smiles,
            "node_label": smiles,
            "uuid": f"reaction_{uuid.uuid4().hex}",
            "yield_predicted": node_data.get("plausibility", 0.0),
            "yield_score": 0,
            "is_predicted": True,
            "is_balanced": False,
            "rxid": smiles,
            "rxsmiles": rxsmiles,
            "node_type": AICP_REACTION_NODE_TYPE,
            "provenance": Provenance(is_in_askcos=True),
            "route_assembly_type": RouteAssemblyType(is_predicted=True, is_evidence=False),
        }

    # Chemical node
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchikey = Chem.MolToInchiKey(mol)
    except Exception:
        raise ValueError(f"RDKit problem parsing SMILES or generating InChIKey: {smiles}")

    return {
        "node_id": smiles,
        "node_label": smiles,
        "uuid": f"substance_{uuid.uuid4().hex}",
        "inchikey": inchikey,
        "canonical_smiles": Chem.CanonSmiles(smiles),
        "srole": "",
        "is_predicted": True,
        "node_type": AICP_SUBSTANCE_NODE_TYPE,
        "provenance": Provenance(is_in_askcos=True),
        "route_assembly_type": RouteAssemblyType(is_predicted=True, is_evidence=False),
    }


def process_uds_edges(
    edges: List[Dict[str, str]],
    node_dict: Dict[str, Any],
    USE_RETRO_RXN_RENDERING: bool = False,
) -> DiGraph:
    """
    Builds a DiGraph from a list of SMILES-keyed edges and a UDS node_dict.

    :param edges: List of {source: smiles, target: smiles} dicts.
    :param node_dict: uds["node_dict"] providing metadata for each SMILES node.
    :param USE_RETRO_RXN_RENDERING: If True, swaps reaction sections for retro rendering.
    :return: DiGraph with AICP-style node and edge metadata.
    """
    graph = DiGraph()
    node_type_map: Dict[str, str] = {}

    # Collect unique node SMILES referenced by edges
    node_smiles: set = set()
    for edge in edges:
        node_smiles.add(edge["source"])
        node_smiles.add(edge["target"])

    # Add nodes
    for smiles in node_smiles:
        if smiles not in node_dict:
            raise ValueError(f"Node SMILES not found in node_dict: {smiles}")
        node_metadata = uds_node_to_synth_node(smiles, node_dict[smiles], USE_RETRO_RXN_RENDERING)
        node_type_map[smiles] = node_metadata["node_type"]
        graph.add_node(smiles, **node_metadata)

    # Add edges (reuses existing convert_askcos_edge_to_synth_edge since edge semantics are identical)
    for edge in edges:
        source, target = edge["source"], edge["target"]
        edge_metadata = convert_askcos_edge_to_synth_edge(source, target, node_type_map, edge_type="reactant_of")
        graph.add_edge(edge_metadata["start_node"], edge_metadata["end_node"], **edge_metadata)

    return graph


def uds_tree2synth_graph(uds: Dict[str, Any], USE_RETRO_RXN_RENDERING: bool = False) -> DiGraph:
    """
    Builds a merged synthesis DiGraph from a UDS dict by combining all pathway edges.
    Equivalent to askcos_tree2synth_paths_with_graph() returning only the graph.

    :param uds: The uds dict from the ASKCOS response (result["uds"]).
    :param USE_RETRO_RXN_RENDERING: If True, swaps reaction sections for retro rendering.
    :return: Merged DiGraph with synthesis roles assigned.
    """
    uuid2smiles: Dict[str, str] = uds["uuid2smiles"]
    pathways: List[List[Dict[str, str]]] = uds["pathways"]
    node_dict: Dict[str, Any] = uds["node_dict"]

    if not pathways:
        raise NoPathsFoundInAskcosResponse("No pathways found in UDS.")

    # Collect unique SMILES-level edges across all pathways
    seen: set = set()
    all_edges: List[Dict[str, str]] = []
    for pathway in pathways:
        for edge in pathway:
            src = uuid2smiles[edge["source"]]
            tgt = uuid2smiles[edge["target"]]
            if (src, tgt) not in seen:
                seen.add((src, tgt))
                all_edges.append({"source": src, "target": tgt})

    synth_graph = process_uds_edges(all_edges, node_dict, USE_RETRO_RXN_RENDERING)
    synth_graph = assign_synth_roles(synth_graph)
    return synth_graph


def uds_tree2synth_paths(uds: Dict[str, Any], USE_RETRO_RXN_RENDERING: bool = False) -> List[Dict[str, Any]]:
    """
    Builds individual synthesis path DiGraphs from a UDS dict.
    Equivalent to askcos_tree2synth_paths() but for the UDS format.

    :param uds: The uds dict from the ASKCOS response (result["uds"]).
    :param USE_RETRO_RXN_RENDERING: If True, swaps reaction sections for retro rendering.
    :return: List of {"idx": int, "path": DiGraph} dicts.
    """
    uuid2smiles: Dict[str, str] = uds["uuid2smiles"]
    pathways: List[List[Dict[str, str]]] = uds["pathways"]
    node_dict: Dict[str, Any] = uds["node_dict"]

    if not pathways:
        raise NoPathsFoundInAskcosResponse("No pathways found in UDS.")

    synth_routes = []
    for idx, pathway in enumerate(pathways):
        edges = [
            {"source": uuid2smiles[e["source"]], "target": uuid2smiles[e["target"]]}
            for e in pathway
        ]
        path_graph = process_uds_edges(edges, node_dict, USE_RETRO_RXN_RENDERING)
        path_graph = assign_synth_roles(path_graph)
        synth_routes.append({"idx": idx, "path": path_graph})

    return synth_routes


def uds_tree2synth_paths_with_graph(
    uds: Dict[str, Any], USE_RETRO_RXN_RENDERING: bool = False
) -> Tuple[DiGraph, List[Dict[str, Any]]]:
    """
    Builds a merged synthesis DiGraph and per-pathway metadata from a UDS dict.
    Equivalent to askcos_tree2synth_paths_with_graph() but for the UDS format.

    :param uds: The uds dict from the ASKCOS response (result["uds"]).
    :param USE_RETRO_RXN_RENDERING: If True, swaps reaction sections for retro rendering.
    :return: Tuple of (merged DiGraph, graph_paths list).
             graph_paths entries: {"path_index": int, "nodes": [node_id, ...], "edges": ["A|B", ...]}
             Edge labels follow the AICP convention of start_node|end_node after direction normalisation,
             which for both UDS edge directions resolves to target|source.
    """
    uuid2smiles: Dict[str, str] = uds["uuid2smiles"]
    pathways: List[List[Dict[str, str]]] = uds["pathways"]
    node_dict: Dict[str, Any] = uds["node_dict"]

    if not pathways:
        raise NoPathsFoundInAskcosResponse("No pathways found in UDS.")

    # Resolve all pathway edges to SMILES up front
    resolved_pathways: List[List[Dict[str, str]]] = []
    for pathway in pathways:
        resolved_pathways.append([
            {"source": uuid2smiles[e["source"]], "target": uuid2smiles[e["target"]]}
            for e in pathway
        ])

    # Build merged graph from all unique SMILES-level edges
    seen: set = set()
    all_edges: List[Dict[str, str]] = []
    for edges in resolved_pathways:
        for edge in edges:
            key = (edge["source"], edge["target"])
            if key not in seen:
                seen.add(key)
                all_edges.append(edge)

    synth_graph = process_uds_edges(all_edges, node_dict, USE_RETRO_RXN_RENDERING)
    synth_graph = assign_synth_roles(synth_graph)

    # Build per-pathway metadata using AICP node/edge identifiers
    # After convert_askcos_edge_to_synth_edge, both edge directions are normalised so that
    # the AICP edge label is always f"{target}|{source}" (start_node|end_node post-flip).
    graph_paths: List[Dict[str, Any]] = []
    for path_index, edges in enumerate(resolved_pathways):
        path_node_ids: set = set()
        path_edge_labels: List[str] = []
        for edge in edges:
            src, tgt = edge["source"], edge["target"]
            path_node_ids.add(src)
            path_node_ids.add(tgt)
            # Edge label after direction normalisation: always tgt|src
            path_edge_labels.append(f"{tgt}|{src}")
        graph_paths.append({
            "path_index": path_index,
            "nodes": list(path_node_ids),
            "edges": path_edge_labels,
        })

    return (synth_graph, graph_paths)
