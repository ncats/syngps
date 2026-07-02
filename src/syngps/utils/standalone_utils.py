from typing import List

from syngps.models import SynthGraph, SynthRoute
from syngps.utils.graph_utils import (
    identify_individual_synthesis_routes,
    grab_yield_predicted,
    project_to_reaction_network,
    annotate_projected_route_with_yield_info,
)
from syngps.yields import aggregate_yields


def find_top_n_routes(synth_graph: SynthGraph, top_n: int) -> List[SynthRoute]:
    """
    Identify and rank the top N synthesis routes by aggregated yield from a SynthGraph.

    Args:
        synth_graph (SynthGraph): A SynthGraph object containing the synthesis graph.
        top_n (int): The maximum number of routes to return.

    Returns:
        List[SynthRoute]: Routes sorted by aggregated_yield descending, capped at top_n.
    """

    viable_routes, _, _ = identify_individual_synthesis_routes(
        synth_graph.synthesis_graph,
        synth_graph.target_molecule_node_id,
    )

    if not viable_routes:
        return []

    predicted_yields = grab_yield_predicted(synth_graph.synthesis_graph)
    for route in viable_routes:
        projected = project_to_reaction_network(route["route_candidate"])
        annotate_projected_route_with_yield_info(projected, predicted_yields)
        route["aggregated_yield"] = aggregate_yields(projected)

    routes_ranked = sorted(
        viable_routes,
        key=lambda r: r.get("aggregated_yield") or 0.0,
        reverse=True,
    )

    return [
        SynthRoute(
            route_candidate=r["route_candidate"],
            route_status=r.get("route_status"),
            method=r.get("method"),
            route_index=r.get("route_index"),
            aggregated_yield=r.get("aggregated_yield"),
        )
        for r in routes_ranked[:top_n]
    ]
