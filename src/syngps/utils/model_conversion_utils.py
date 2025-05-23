from networkx import DiGraph
from syngps.models import (
    TopNYieldSynthRoutesResult,
    SynthGraph,
    SynthRoute,
    SynthView,
    SynthRouteContainer,
    CombinationGraph,
    SynthesisRoutesResponse
)
from syngps.utils.graph_json_utils import json_to_graph, create_subgraphs_from_routes
from typing import Any


def build_synthgraph_from_json(synth_graph_json: dict) -> SynthGraph:
    # Extract graph and convert to networkx
    synthesis_graph: DiGraph = json_to_graph(synth_graph_json)

    # Extract fields for the model
    search_params = synth_graph_json.get("search_params")
    time_info = synth_graph_json.get("time_info")
    summary = synth_graph_json.get("summary")
    availability = synth_graph_json.get("availability")  # Optional

    # Infer the target molecule node ID (srole = 'tm')
    target_molecule_node_id = _extract_target_node_label(synthesis_graph)

    # Raise an error if the target molecule node ID is not found
    if not target_molecule_node_id:
        raise ValueError("Target molecule node ID not found in synthesis graph")

    # Create the SynthGraph instance
    sg = SynthGraph(
        synthesis_graph=synthesis_graph,
        search_params=search_params,
        time_info=time_info,
        summary=summary,
        availability=availability,
        target_molecule_node_id=target_molecule_node_id,
    )
    return sg


def build_topn_yield_result(json_dict: dict) -> TopNYieldSynthRoutesResult:
    # Step 1: Extract and convert synthesis_graph
    synth_graph_json = json_dict.get("synth_graph")

    if not synth_graph_json:
        raise ValueError("Synthesis graph not found in JSON data")

    synthesis_graph = json_to_graph(synth_graph_json)

    # Raise an error if synthesis_graph is not found, it is needed to reconstruct subgraphs
    if synthesis_graph is None:
        raise ValueError("Unable to extract synthesis graph from JSON data")

    # Extract target molecule node ID
    target_node_id = _extract_target_node_label(synthesis_graph)
    if not target_node_id:
        raise ValueError("Target molecule node ID not found in synthesis graph")

    # Step 2: Create SynthGraph instance
    synth_graph = SynthGraph(
        synthesis_graph=synthesis_graph,
        search_params=json_dict["search_params"],
        target_molecule_node_id=target_node_id,
        time_info=json_dict.get("time_info"),
        summary=json_dict.get("summary"),
        availability=json_dict.get("availability"),
    )

    # Step 3: Build artifact_free_graph
    artifact_free_graph = SynthView(view=synthesis_graph)

    # Step 4: Build list of SynthRoutes
    if "routes" not in json_dict:
        raise ValueError("Routes not found in JSON data")
    routes = create_subgraphs_from_routes(json_dict["routes"], synthesis_graph)

    # Step 5: Build route_candidates (if present)
    route_candidates = None
    if "route_candidates" in json_dict and json_dict["route_candidates"]:
        subgraphs = []
        for route_data in json_dict["route_candidates"].get("subgraphs", []):
            route_graph = json_to_graph(route_data)
            subgraphs.append(SynthRoute(
                route_candidate=route_graph,
                aggregated_yield=route_data.get("aggregated_yield"),
                predicted=route_data.get("predicted", False),
                source=route_data.get("source"),
                route_index=route_data.get("route_index"),
                route_status=route_data.get("route_status"),
                method=route_data.get("method"),
            ))
        route_candidates = SynthRouteContainer(subgraphs=subgraphs)

    # Step 6: Build combination_graphs (if present)
    combination_graphs = None
    if "combination_graphs" in json_dict:
        combination_graphs = []
        for cg_data in json_dict["combination_graphs"]:
            cg_graph = json_to_graph(cg_data)
            combination_graphs.append(CombinationGraph(
                route_index=cg_data.get("route_index"),
                route_candidate=cg_graph
            ))

    # Step 7: Build the final TopNYieldSynthRoutesResult
    result = TopNYieldSynthRoutesResult(
        search_params=json_dict["search_params"],
        synth_graph=synth_graph,
        artifact_free_graph=artifact_free_graph,
        combination_graphs=combination_graphs,
        route_candidates=route_candidates,
        routes=routes,
        time_info=json_dict["time_info"],
        summary=json_dict["summary"],
        availability=json_dict.get("availability"),
    )

    return result


def build_synth_routes_result(json_dict: dict) -> SynthesisRoutesResponse:
    # Step 1: Extract evidence graph
    if "evidence_synth_graph" in json_dict:
        evidence_synth_graph = build_synthgraph_from_json(json_dict["evidence_synth_graph"])
    else:
        evidence_synth_graph = None

    # Step 2: Extract predictive graph
    if "predictive_synth_graph" in json_dict:
        predictive_synth_graph = build_synthgraph_from_json(json_dict["predictive_synth_graph"])  
    else:
        predictive_synth_graph = None

    # Step 3: Validate at least one graph is present
    if evidence_synth_graph is None and predictive_synth_graph is None:
        raise ValueError("Neither evidence_synth_graph nor predictive_synth_graph found in JSON data")
    
    # Step 4: Create evidence routes
    evidence_routes_mapped = []
    if evidence_synth_graph:
        evidence_routes = [route for route in json_dict["routes"] if not route.get("predicted", False)]
        evidence_routes_mapped = create_subgraphs_from_routes(evidence_routes, evidence_synth_graph.synthesis_graph)

    # Step 5: Create predictive routes
    predictive_routes_mapped = []
    if predictive_synth_graph:
        predictive_routes = [route for route in json_dict["routes"] if route.get("predicted", False)]
        predictive_routes_mapped = create_subgraphs_from_routes(predictive_routes, predictive_synth_graph.synthesis_graph)

    # Step 6: Build the final SynthesisRoutesResponse
    result = SynthesisRoutesResponse(
        target_molecule_inchikey=json_dict["target_molecule_inchikey"],
        target_molecule_smiles=json_dict["target_molecule_smiles"],
        reaction_steps=json_dict["reaction_steps"],
        evidence_routes_success=json_dict["evidence_routes_success"],
        predicted_routes_success=json_dict["predicted_routes_success"],
        routes=evidence_routes_mapped + predictive_routes_mapped,
        evidence_synth_graph=evidence_synth_graph,
        predictive_synth_graph=predictive_synth_graph,
    )

    return result


def _extract_target_node_label(G: DiGraph) -> str | None:
    """
    Retrieves the label of the target molecule node in the synthesis graph.
    """
    for node, data in G.nodes(data=True):
        if data.get("srole") == "tm":
            return node
    return None
