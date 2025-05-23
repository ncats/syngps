import cProfile
import pstats
import time
from typing import Any, Dict, List, Tuple

from networkx import DiGraph
import yields as ay
from adapters.abstract_data_adapter import AbstractDataAdapter
from errors.errors import (
    SynthGraphParsingException,
)
from models.models import (
    ProjectedRoute,
    SGPInput,
    SynthGraph,
    SynthGraphSearch,
    SynthRoute,
    SynthRouteContainer,
    SynthView,
    TimeInfo,
    TopNYieldSearch,
    TopNYieldSynthRoutesResult,
    TopologyInfo,
)
from utils import graph_json_utils, graph_utils
import logging
logger = logging.getLogger(__name__)


class AicpFunctions:
    """
    Core AICP functionality. This class is responsible for fetching synthesis graph and finding top N yield synthesis routes.
    This class gets initialized with data adapter. It uses these adapters to fetch data from the knowledge base.
    """

    def __init__(
        self,
        data_adapter: AbstractDataAdapter,
        enable_profiler: bool = False,
    ):
        """
        Initializes AicpFunctions with data adapter and inventory adapter.

        Args:
            data_adapter (AbstractDataAdapter): Data adapter to fetch data from the knowledge base.
            enable_profiler (bool, optional): Enable profiler. Defaults to False.
        """
        self.enable_profiler = enable_profiler
        self.profiler = cProfile.Profile()
        self.data_adapter: AbstractDataAdapter = data_adapter


    def fetch_synthesis_graph(self, synth_search: SynthGraphSearch) -> SynthGraph:
        """
        Fetches synthesis graph for the given search parameters.
        Core functionality is found in _fetch_synthesis_graph method, which is called by this method with a profiler wrapper.

        Args:
            synth_search (SynthGraphSearch): Search parameters to fetch synthesis graph.

        Returns:
            SynthGraph: Synthesis graph for the given search parameters.
        """
        if self.enable_profiler:
            logger.debug("Enabling profiler for fetch_synthesis_graph...")
            self.profiler.enable()

        synth_graph, synth_graph_time, artifact_free_graph_obj = self._fetch_synthesis_graph(synth_search)

        if self.enable_profiler:
            logger.debug("Stopping profiler for fetch_synthesis_graph...")
            self.profiler.disable()
            ps = pstats.Stats(self.profiler).sort_stats("cumulative")
            synth_graph.profiling_stats = ps
            logger.debug(f"Profiling stats:\n{ps}")

        return synth_graph, artifact_free_graph_obj

    def find_top_n_yield_synthesis_routes(self, route_search: TopNYieldSearch, top_n: int) -> TopNYieldSynthRoutesResult:
        """
        Finds top N yield synthesis routes for the given search parameters.
        Core functionality is found in _find_top_n_yield_synthesis_routes method, which is called by this method with a profiler wrapper.

        Args:
            route_search (TopNYieldSearch): Search parameters to find top N yield synthesis routes.
            top_n (int): Number of top yield synthesis routes to find.

        Returns:
            TopNYieldSynthRoutesResult: Top N yield synthesis routes for the given search parameters.
        """
        if self.enable_profiler:
            logger.debug("Enabling profiler for find_top_n_yield_synthesis_routes...")
            self.profiler.enable()

        top_n_yield_synthesis_routes = self._find_top_n_yield_synthesis_routes(route_search, top_n)

        if self.enable_profiler:
            logger.debug("Stopping profiler for find_top_n_yield_synthesis_routes...")
            self.profiler.disable()
            ps = pstats.Stats(self.profiler).sort_stats("cumulative")
            top_n_yield_synthesis_routes.profiling_stats = ps
            logger.debug(f"Profiling stats:\n{ps}")

        return top_n_yield_synthesis_routes

    def prune_synth_graph(self, sgp_input: SGPInput) -> Dict:
        """
        Prunes the synthesis graph based on the search parameters.
        """

        try:
            synth_graph = graph_json_utils.json_to_graph(sgp_input.synth_graph_json.model_dump())
        except Exception:
            logger.error("Error parsing synthesis graph JSON:", exc_info=True)
            raise SynthGraphParsingException("Error parsing synthesis graph JSON")

        unwanted_substances = sgp_input.unwanted_substances
        unwanted_reactions = sgp_input.unwanted_reactions

        logger.debug(f"Pruning graph with: unwanted_substances={unwanted_substances}, unwanted_reactions={unwanted_reactions}")
        nodes_to_remove = []
        # Mark all nodes to remove
        for node, attrs in synth_graph.nodes(data=True):
            if attrs.get("node_type", "").lower() == "substance":
                if attrs.get("inchikey") in unwanted_substances:
                    nodes_to_remove.append(attrs.get("inchikey"))
            elif attrs.get("node_type", "").lower() == "reaction":
                if attrs.get("rxid") in unwanted_reactions:
                    nodes_to_remove.append(attrs.get("rxid"))

        # Remove nodes
        synth_graph.remove_nodes_from(nodes_to_remove)
        logger.debug(f"Removed {len(nodes_to_remove)} nodes from the graph.")

        # Remove orphan nodes
        orphan_nodes = [node for node in synth_graph.nodes() if synth_graph.degree(node) == 0]
        synth_graph.remove_nodes_from(orphan_nodes)
        logger.debug(f"Removed {len(orphan_nodes)} orphan nodes from the graph.")

        logger.debug(f"Final Graph size: {len(synth_graph.nodes)} nodes, {len(synth_graph.edges)} edges.")
        return graph_json_utils.graph_to_json(synth_graph)

    def _fetch_synthesis_graph(self, route_search: TopNYieldSearch, synthesis_graph_json=None) -> SynthGraph:
        """
        Fetches synthesis graph for the given search parameters.

        Args:
            route_search (TopNYieldSearch): Search parameters to find top N yield synthesis routes.
            synthesis_graph_json: Optionally provide the synthesis_graph_json and skip fetch synthesis graph

        Returns:
            SynthGraph: Synthesis graph for the given search parameters.
        """
        synth_graph_time_start = time.time()
        if synthesis_graph_json:
            # Manually construct a directed graph (DiGraph) from JSON
            G = DiGraph()

            # Add nodes
            for node in route_search.synthesis_graph_json.get("nodes", []):
                G.add_node(node["node_label"], **node)

            # Add edges
            for edge in route_search.synthesis_graph_json.get("edges", []):
                G.add_edge(edge["start_node"], edge["end_node"], **edge)

            target_node = route_search.target_molecule_inchikey
            synth_graph = SynthGraph(
                target_molecule_node_id=target_node,
                synthesis_graph=G.copy(),
                search_params=route_search,
            )
        else:
            synth_graph = self.data_adapter.fetch_synthesis_graph(route_search)

            # Assign inventory roles to substances
            self._assign_substance_inventory_roles(synth_graph, route_search.target_molecule_inchikey)

        # Remove incompatible reactions
        artifact_free_synth_graph = graph_utils.remove_incompatible_reactions(
            synth_graph.synthesis_graph,
            synth_graph.target_molecule_node_id,
        )
        synth_graph_time = time.time() - synth_graph_time_start
        artifact_free_graph_obj = SynthView(view=artifact_free_synth_graph)

        return synth_graph, synth_graph_time, artifact_free_graph_obj

    def _find_top_n_yield_synthesis_routes(self, route_search: TopNYieldSearch, top_n: int) -> TopNYieldSynthRoutesResult:
        """
        Finds top N yield synthesis routes for the given search parameters.

        Args:
            route_search (TopNYieldSearch): Search parameters to find top N yield synthesis routes.
            top_n (int): Number of top yield synthesis routes to find.

        Returns:
            TopNYieldSynthRoutesResult: Top N yield synthesis routes for the given search parameters.
        """

        # Create a synthesis graph from search params
        synthesis_graph_json = route_search.synthesis_graph_json
        synth_graph, synth_graph_time, artifact_free_graph = self._fetch_synthesis_graph(route_search, synthesis_graph_json=synthesis_graph_json)

        enumerate_routes_time_start = time.time()

        # Identify synthesis routes
        routes, route_candidates, combination_graphs = self._identify_synth_graph_routes(
            synth_graph, route_search.include_route_candidates, route_search.include_combination_graphs
        )

        enumerate_routes_time = time.time() - enumerate_routes_time_start
        topology_info = TopologyInfo(num_routes=len(routes), num_combination_graphs=len(combination_graphs), num_route_candidates=len(route_candidates))
        # TODO : Find better way to validate the graph
        if artifact_free_graph.view is None:
            raise ValueError("'synth_graph' is missing synthesis_graph property.")

        G_dag = artifact_free_graph.view.copy()

        individual_synthesis_routes: List[SynthRoute] = routes

        aggregate_yield_time_start = time.time()

        # Get the predicted yields for the synthesis graph
        yields = graph_utils.grab_yield_predicted(G_dag)
        projected_routes: List[ProjectedRoute] = []
        # Project the synthesis graph to a reaction network and annotate it with the predicted yields
        for route in individual_synthesis_routes:
            idx = route.route_index
            # Project the synthesis tree to a reaction network
            projected_route: DiGraph = graph_utils.project_to_reaction_network(route.route_candidate)
            # Annotate the projected route with the predicted yields
            graph_utils.annotate_projected_route_with_yield_info(projected_route, yields)
            # Aggregate the yields of the projected route
            route_aggregated_yield = ay.aggregate_yields(projected_route)
            # Add the projected route to the list of projected routes
            projected_routes.append(
                ProjectedRoute(
                    route_index=idx,
                    route=projected_route,
                    aggregated_yield=route_aggregated_yield,
                )
            )

        aggregate_yield_time = time.time() - aggregate_yield_time_start

        rank_routes_time_start = time.time()

        final_projected_routes: List[ProjectedRoute] = self._top_n_yield_routes(projected_routes, top_n)

        projected_route_indices = {pr.route_index for pr in final_projected_routes}

        # Filter out routes in synth_graph.routes that don't correspond to the selected projected routes
        routes = [route for route in routes if route.route_index in projected_route_indices]
        rank_routes_time = time.time() - rank_routes_time_start

        tree_combination: Dict[int, Any] = {}
        for proute in final_projected_routes:
            tree_combination[proute.route_index] = {"projected_route": proute}

        for sroute in individual_synthesis_routes:
            if sroute.route_index in tree_combination.keys():
                tree_combination[sroute.route_index]["synth_route"] = sroute

        for idx, tree in tree_combination.items():
            pr: ProjectedRoute = tree["projected_route"]
            aggregated_yield = pr.aggregated_yield if pr.aggregated_yield is not None else 0.0

            # Find the corresponding route in synth_graph.routes based on the route_index
            for route in routes:
                if route.route_index == pr.route_index:
                    route.aggregate_yield = aggregated_yield
                    break  # Once you find the route, break out of the loop

        # Create the final TopNYieldSynthRoutesResult object
        top_n_results: TopNYieldSynthRoutesResult = TopNYieldSynthRoutesResult(
            search_params=route_search,
            synth_graph=synth_graph,
            artifact_free_graph=artifact_free_graph,
            combination_graphs=SynthRouteContainer(subgraphs=combination_graphs),
            route_candidates=SynthRouteContainer(subgraphs=route_candidates),
            routes=SynthRouteContainer(subgraphs=routes),
            time_info=TimeInfo(
                synth_graph_time=synth_graph_time,
                enumerate_routes_time=enumerate_routes_time,
                aggregate_yield_time=aggregate_yield_time,
                rank_routes_time=rank_routes_time,
            ),
            summary=topology_info,
        )
        return top_n_results

    def _assign_substance_inventory_roles(self, synth_graph: SynthGraph, target_molecule_inchikey: str) -> None:
        """
        Assigns inventory roles to substances in the given synthesis graph using the inventory adapter.

        Args:
            synth_graph (SynthGraph): Synthesis graph to assign inventory roles to substances.
        """
        substances = [node_id for node_id, node_data in synth_graph.synthesis_graph.nodes(data=True) if node_data.get("node_type") == "substance"]

        for substance in substances:
            inchikey = synth_graph.synthesis_graph.nodes[substance]["inchikey"]
            in_degree = synth_graph.synthesis_graph.in_degree(substance)
            if inchikey == target_molecule_inchikey:
                srole = "tm"
            # elif self.inventory_adapter.in_inventory_by_inchikey(inchikey):
            #     substance.is_in_inventory = True
            #     srole = "sm"
            elif in_degree == 0 and synth_graph.search_params.leaves_as_sm:
                srole = "sm"
            else:
                srole = "im"
            synth_graph.synthesis_graph.nodes[substance]["srole"] = srole

    def _identify_synth_graph_routes(
        self, synth_graph: SynthGraph, include_route_candidates: bool, include_combination_graphs: bool
    ) -> Tuple[List[SynthRoute], List[SynthRoute], List[SynthRoute]]:
        """
        Identifies synthesis routes for the given synthesis graph and adds them to the SynthGraph object.

        Args:
            synth_graph (SynthGraph): Synthesis graph to identify synthesis routes for.
        """
        route_candidates: List[SynthRoute] = []

        individual_synthesis_routes, individual_synthesis_routes_candidates, combination_graphs = graph_utils.identify_individual_synthesis_routes(
            synth_graph.synthesis_graph,
            synth_graph.target_molecule_node_id,
        )

        routes = [SynthRoute(**route) for route in individual_synthesis_routes]
        if include_route_candidates:
            route_candidates = [SynthRoute(**route) for route in individual_synthesis_routes_candidates]
        if include_combination_graphs:
            combination_maps = []
            for idx, combination_graph in enumerate(combination_graphs, start=1):
                cg_map = {}
                cg_map["route_candidate"] = combination_graph
                cg_map["route_index"] = idx
                combination_maps.append(cg_map)
            combination_graphs_parsed = [SynthRoute(**route) for route in combination_maps]
        else:
            combination_graphs_parsed = []
        return routes, route_candidates, combination_graphs_parsed

    def _top_n_yield_routes(self, all_routes_with_aggregated_yields: List[ProjectedRoute], top_n: int) -> List[ProjectedRoute]:
        """
        Returns the top N routes with the highest aggregated yield values.

        This function sorts the given list of ProjectedRoute objects in descending order based on their aggregated yields
        and returns the top N routes. Routes with a None value for aggregated_yield are considered the lowest and appear last.

        Args:
            all_routes_with_aggregated_yields (List[ProjectedRoute]): A list of ProjectedRoute objects to be sorted.
            top_n (int): The number of top yielding routes to return. Must be greater than 0.

        Returns:
            List[ProjectedRoute]: A list containing the top N routes sorted by their aggregated yields.

        Raises:
            ValueError: If top_n is less than 0.
        """
        if top_n < 1:
            raise ValueError("top_n must be greater than 0.")

        sorted_routes = sorted(
            all_routes_with_aggregated_yields,
            key=lambda x: x.aggregated_yield if x.aggregated_yield is not None else 0.0,
            reverse=True,
        )

        return sorted_routes[:top_n]

    def aggregate_route_yield(self, synth_route: SynthRoute) -> SynthRoute:
        """
        Aggregates the yield of a given synthesis route
        """
        # Get the predicted yields for the synthesis graph
        yields = graph_utils.grab_yield_predicted(synth_route.route_candidate)

        # Project route to a reaction network
        projected_route: DiGraph = graph_utils.project_to_reaction_network(synth_route.route_candidate)
        # Ensure projected route has yield information
        graph_utils.annotate_projected_route_with_yield_info(projected_route, yields)
        # Calculate ag yield
        route_aggregated_yield = ay.aggregate_yields(projected_route)
        # Return route with ag yield
        synth_route.aggregate_yield = route_aggregated_yield
        return synth_route
