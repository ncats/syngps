import requests
import networkx as nx
import logging
from syngps.models import (
    SynthGraph,
    TopNYieldSynthRoutesResult,
    SynthRoute,
    SynthesisRoutesResponse
)
import uuid
from typing import Optional


# Setup logger
logger = logging.getLogger(__name__)

# Constants
CYTOSCAPE_API_URL = "http://localhost:1234/v1"
LAYOUT_TYPE = "hierarchical"
STYLE_NAME = "New SynGPS API"

# Custom Exceptions


class CytoscapeUnavailableError(Exception):
    pass


class NetworkCreationError(Exception):
    pass


class NetworkViewError(Exception):
    pass


class StyleApplicationError(Exception):
    pass


class LayoutApplicationError(Exception):
    pass


STYLE_JSON = {
    "title": STYLE_NAME,
    "defaults": [
        {"visualProperty": "NODE_SIZE", "value": 40},
        {"visualProperty": "EDGE_LINE_TYPE", "value": "SOLID"},
        {"visualProperty": "EDGE_WIDTH", "value": 2},
        {"visualProperty": "EDGE_CURVED", "value": False},
        {"visualProperty": "EDGE_TARGET_ARROW_SHAPE", "value": "DELTA"},
    ],
    "mappings": [
        {
            "mappingType": "discrete",
            "mappingColumn": "srole",
            "mappingColumnType": "String",
            "visualProperty": "NODE_FILL_COLOR",
            "map": [
                {"key": "tm", "value": "#4C8DA6"},
                {"key": "im", "value": "#AAAAAA"},
                {"key": "sm", "value": "#D8C571"},
            ],
        },
        {
            "mappingType": "discrete",
            "mappingColumn": "node_type",
            "mappingColumnType": "String",
            "visualProperty": "NODE_SHAPE",
            "map": [
                {"key": "substance", "value": "ROUND_RECTANGLE"},
                {"key": "reaction", "value": "ELLIPSE"},
            ],
        },
        {
            "mappingType": "discrete",
            "mappingColumn": "edge_type",
            "mappingColumnType": "String",
            "visualProperty": "EDGE_STROKE_UNSELECTED_PAINT",
            "map": [
                {"key": "product_of", "value": "#EC7014"},
                {"key": "reactant_of", "value": "#225EA8"},
                {"key": "reagent_of", "value": "#225EA8"},
            ],
        },
        {
            "mappingType": "discrete",
            "mappingColumn": "edge_type",
            "mappingColumnType": "String",
            "visualProperty": "EDGE_TARGET_ARROW_UNSELECTED_PAINT",
            "map": [
                {"key": "product_of", "value": "#EC7014"},
                {"key": "reactant_of", "value": "#225EA8"},
                {"key": "reagent_of", "value": "#225EA8"},
            ],
        },
        {
            "mappingType": "passthrough",
            "mappingColumn": "node_label",
            "mappingColumnType": "String",
            "visualProperty": "NODE_LABEL",
        },
    ],
}


def validate_cytoscape(timeout: float = 1.0) -> bool:
    """Returns True if Cytoscape REST API is available, else raises exception."""
    try:
        response = requests.get(f"{CYTOSCAPE_API_URL}/version", timeout=timeout)
        if not response.ok:
            raise CytoscapeUnavailableError("Cytoscape REST API responded with error.")
        logger.debug("Cytoscape REST API is available.")
        return True
    except requests.RequestException as e:
        raise CytoscapeUnavailableError(f"Cannot connect to Cytoscape REST API: {e}")


def to_cytoscape_json(G: nx.DiGraph) -> dict:
    """Converts a DiGraph to Cytoscape JSON format."""
    # Step 1: Collect nodes to rename
    nodes_to_rename = []
    for n in G.nodes:
        new_id = G.nodes[n].get('uuid')
        G.nodes[n]['node_id'] = new_id
        if new_id is not None and n != new_id:
            nodes_to_rename.append((n, new_id))

    # Step 2: Perform renaming
    for old_id, new_id in nodes_to_rename:
        rename_node(G, old_id, new_id)

    # Step 3: Convert to Cytoscape JSON
    return nx.cytoscape_data(G, name="node_label", ident="node_id")


def rename_node(G: nx.DiGraph, old_id, new_id):
    if old_id not in G:
        raise ValueError(f"Node {old_id} not found in graph.")
    if new_id in G:
        raise ValueError(f"Node {new_id} already exists in graph.")

    # Add the new node with the same attributes
    G.add_node(new_id, **G.nodes[old_id])

    # Redirect incoming edges
    for u, _, data in G.in_edges(old_id, data=True):
        G.add_edge(u, new_id, **data)

    # Redirect outgoing edges
    for _, v, data in G.out_edges(old_id, data=True):
        G.add_edge(new_id, v, **data)

    # Remove the old node
    G.remove_node(old_id)


def send_to_cytoscape(cyjs: dict) -> int:
    """Uploads the JSON network to Cytoscape and returns its SUID."""
    validate_cytoscape()

    resp = requests.post(f"{CYTOSCAPE_API_URL}/networks?format=cyjs", json=cyjs)
    if not resp.ok:
        logger.error("Failed to create network.")
        raise NetworkCreationError(resp.text)

    suid = resp.json().get("networkSUID")
    logger.debug(f"Network created with SUID: {suid}")

    view_resp = requests.get(f"{CYTOSCAPE_API_URL}/networks/{suid}/views/first")
    if not view_resp.ok:
        logger.error("Failed to create network view.")
        raise NetworkViewError(view_resp.text)

    logger.debug("Network view created.")
    return int(view_resp.json()["data"]["SUID"])


def ensure_style() -> None:
    """Ensures the custom style exists in Cytoscape."""
    validate_cytoscape()

    resp = requests.get(f"{CYTOSCAPE_API_URL}/styles")
    if not resp.ok:
        raise StyleApplicationError("Failed to retrieve styles.")

    if STYLE_NAME in resp.json():
        logger.debug(f"Style '{STYLE_NAME}' already exists.")
        return

    logger.debug(f"Creating style '{STYLE_NAME}'...")
    create_resp = requests.post(f"{CYTOSCAPE_API_URL}/styles", json=STYLE_JSON)
    if not create_resp.ok:
        raise StyleApplicationError("Failed to create new style.")

    logger.debug(f"Style '{STYLE_NAME}' created.")


def apply_style(suid: int) -> None:
    """Applies the custom style to the given network."""
    validate_cytoscape()
    ensure_style()

    resp = requests.get(f"{CYTOSCAPE_API_URL}/apply/styles/{STYLE_NAME}/{suid}")
    if not resp.ok:
        logger.error(f"Failed to apply style '{STYLE_NAME}'.")
        raise StyleApplicationError(resp.text)

    logger.debug(f"Style '{STYLE_NAME}' applied to network {suid}.")


def apply_layout(suid: int) -> None:
    """Applies the layout to the given network."""
    validate_cytoscape()

    resp = requests.get(f"{CYTOSCAPE_API_URL}/apply/layouts/{LAYOUT_TYPE}/{suid}")
    if not resp.ok:
        logger.error(f"Failed to apply layout '{LAYOUT_TYPE}'.")
        raise LayoutApplicationError(resp.text)

    logger.debug(f"Layout '{LAYOUT_TYPE}' applied to network {suid}.")


def send_synthgraph_to_cytoscape(sg: SynthGraph, collection_name: Optional[str] = None) -> int:
    """
    Sends a SynthGraph to Cytoscape, assigns it to a collection,
    and applies style and layout.

    Args:
        sg (SynthGraph): The synthesis graph to send.
        collection_name (str, optional): Name of the collection to group the network into.

    Returns:
        int: The SUID of the created network.
    """
    logger.info("Sending SynthGraph to Cytoscape...")

    # Convert DiGraph to Cytoscape JSON
    cyjs = to_cytoscape_json(sg.synthesis_graph)

    # Assign default collection name if not provided
    if not collection_name:
        target = sg.target_molecule_node_id or "synth"
        collection_name = f"SynthGraph_{target}"

    # Inject network-level metadata
    cyjs["data"] = {
        "name": collection_name,
        "networkCollection": collection_name
    }

    # Send to Cytoscape
    suid = send_to_cytoscape(cyjs)

    # Apply style and layout
    apply_style(suid)
    apply_layout(suid)

    logger.info(f"SynthGraph sent and styled in Cytoscape with SUID {suid}.")
    return suid


def send_network_to_cytoscape(network: nx.DiGraph, network_name: Optional[str] = None) -> int:
    """
    Sends a network to Cytoscape, assigns it to a collection,
    and applies style and layout.

    Args:
        network (nx.DiGraph): The network to send.
        network_name (str, optional): Name of the network to group the network into.

    Returns:
        int: The SUID of the created network.
    """
    logger.info("Sending Network to Cytoscape...")

    # Convert DiGraph to Cytoscape JSON
    cyjs = to_cytoscape_json(network)

    # Assign default collection name if not provided
    if not network_name:
        network_name = "Network"

    # Inject network-level metadata
    cyjs["data"] = {
        "name": network_name,
        "networkCollection": network_name
    }

    # Send to Cytoscape
    suid = send_to_cytoscape(cyjs)

    # Apply style and layout
    apply_style(suid)
    apply_layout(suid)

    logger.info(f"Network sent and styled in Cytoscape with SUID {suid}.")
    return suid

def send_synthgraphs_to_cytoscape(synthgraphs: list[SynthGraph], collection_name: Optional[str] = None) -> list[int]:
    """
    Sends a list of SynthGraphs to Cytoscape.
    """
    return [send_synthgraph_to_cytoscape(sg, collection_name) for sg in synthgraphs]


def send_synthroute_to_cytoscape(synthroute: SynthRoute, collection_name: Optional[str] = None) -> int:
    """
    Sends a SynthRoute to Cytoscape, assigns it to a collection,
    and applies style and layout.

    Args:
        synthroute (SynthRoute): The synthesis graph to send.
        collection_name (str, optional): Name of the collection to group the network into.

    Returns:
        int: The SUID of the created network.
    """
    logger.info("Sending SynthRoute to Cytoscape...")

    # Convert DiGraph to Cytoscape JSON
    cyjs = to_cytoscape_json(synthroute.route_candidate)

    # Assign default collection name if not provided
    if not collection_name:
        collection_name = f"SynthRoute_{"Predicted" if synthroute.predicted else "Evidence"}_{synthroute.route_index}"

    # Inject network-level metadata
    cyjs["data"] = {
        "name": collection_name,
        "networkCollection": collection_name,
        "aggregated_yield": synthroute.aggregated_yield if synthroute.aggregated_yield else "N/A"
    }

    # Send to Cytoscape
    suid = send_to_cytoscape(cyjs)

    # Apply style and layout
    apply_style(suid)
    apply_layout(suid)

    logger.info(f"SynthRoute sent and styled in Cytoscape with SUID {suid}.")
    return suid


def send_synthroutes_to_cytoscape(synthroutes: list[SynthRoute], collection_name: Optional[str] = None) -> list[int]:
    """
    Sends a list of SynthRoutes to Cytoscape.
    """
    return [send_synthroute_to_cytoscape(sr, collection_name) for sr in synthroutes]


def send_topn_results_to_cytoscape(topn_results: TopNYieldSynthRoutesResult, collection_name: Optional[str] = None) -> list[int]:
    """
    Sends a a top N yield synthesis routes result to Cytoscape.
    """

    # Send synth graph
    sg_suid = send_synthgraph_to_cytoscape(topn_results.synth_graph, collection_name)

    # Send synth routes
    sr_suids = send_synthroutes_to_cytoscape(topn_results.routes, collection_name)

    return [sg_suid] + sr_suids


def send_synthesis_routes_response_to_cytoscape(synthesis_routes_response: SynthesisRoutesResponse, collection_name: Optional[str] = None) -> list[int]:
    """
    Sends a SynthesisRoutesResponse to Cytoscape.
    """

    sg_suids = []

    # Send evidence synth graph
    if synthesis_routes_response.evidence_synth_graph:
        esg_suid = send_synthgraph_to_cytoscape(synthesis_routes_response.evidence_synth_graph, collection_name)
        sg_suids.append(esg_suid)

    # Send predictive synth graph
    if synthesis_routes_response.predictive_synth_graph:
        psg_suid = send_synthgraph_to_cytoscape(synthesis_routes_response.predictive_synth_graph, collection_name)
        sg_suids.append(psg_suid)

    # Filter out routes that are not SynthRoute objects
    routes_to_send = []
    for route in synthesis_routes_response.routes:
        if isinstance(route, SynthRoute):
            routes_to_send.append(route)
        else:
            logger.info(
                f"Skipping non-SynthRoute object of type {type(route).__name__} when sending to Cytoscape."
            )

    # Send synth routes
    sr_suids = send_synthroutes_to_cytoscape(routes_to_send, collection_name)

    return sg_suids + sr_suids
