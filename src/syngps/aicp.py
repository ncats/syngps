import cProfile
import pstats
import time
from typing import Any, Dict, List, Tuple
from networkx import DiGraph
from syngps import yields as ay
from syngps.adapters.abstract_data_adapter import AbstractDataAdapter
from syngps.errors import (
    SynthGraphParsingException,
)
from syngps.models.models import (
    Availability,
    CombinationGraph,
    CommercialAvailability,
    Inventory,
    ProjectedRoute,
    SGPInput,
    SynthGraph,
    SynthGraphJson,
    SynthGraphSearch,
    SynthRoute,
    SynthRouteContainer,
    SynthView,
    TimeInfo,
    TopNYieldSearch,
    TopNYieldSynthRoutesResult,
    TopologyInfo,
)
from syngps.utils import graph_json_utils, graph_utils
import logging

# Setup logger
logger = logging.getLogger(__name__)


class AicpFunctions:
    """
    Core AICP functionality. This class is responsible for fetching synthesis graph and finding top N yield synthesis routes.
    This class gets initialized with data adapter and inventory adapter. It uses these adapters to fetch data from the knowledge base and inventory.
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

        synth_graph = self._fetch_synthesis_graph(synth_search)

        if self.enable_profiler:
            logger.debug("Stopping profiler for fetch_synthesis_graph...")
            self.profiler.disable()
            ps = pstats.Stats(self.profiler).sort_stats("cumulative")
            synth_graph.profiling_stats = ps
            logger.debug(f"Profiling stats:\n{ps}")

        return synth_graph

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
            logger.debug(
                "Enabling profiler for find_top_n_yield_synthesis_routes...")
            self.profiler.enable()

        top_n_yield_synthesis_routes = self._find_top_n_yield_synthesis_routes(
            route_search, top_n)

        if self.enable_profiler:
            logger.debug(
                "Stopping profiler for find_top_n_yield_synthesis_routes...")
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
            synth_graph = graph_json_utils.json_to_graph(
                sgp_input.synth_graph_json.model_dump())
        except Exception:
            logger.error("Error parsing synthesis graph JSON:", exc_info=True)
            raise SynthGraphParsingException(
                "Error parsing synthesis graph JSON")

        unwanted_substances = sgp_input.unwanted_substances
        unwanted_reactions = sgp_input.unwanted_reactions

        logger.debug(
            f"Pruning graph with: unwanted_substances={unwanted_substances}, unwanted_reactions={unwanted_reactions}")
        nodes_to_remove = []
        # Mark all nodes to remove
        for node, attrs in synth_graph.nodes(data=True):
            if attrs.get("node_type", "").lower() == "substance":
                logger.info("Node type: Substance")
                logger.info(attrs)
                if attrs.get("inchikey") in unwanted_substances:
                    nodes_to_remove.append(attrs.get("inchikey"))
            elif attrs.get("node_type", "").lower() == "reaction":
                if attrs.get("rxid") in unwanted_reactions:
                    nodes_to_remove.append(attrs.get("rxid"))

        # Remove nodes
        synth_graph.remove_nodes_from(nodes_to_remove)
        logger.debug(f"Removed {len(nodes_to_remove)} nodes from the graph.")

        # Remove orphan nodes
        orphan_nodes = [node for node in synth_graph.nodes()
                        if synth_graph.degree(node) == 0]
        synth_graph.remove_nodes_from(orphan_nodes)
        logger.debug(
            f"Removed {len(orphan_nodes)} orphan nodes from the graph.")

        logger.debug(
            f"Final Graph size: {len(synth_graph.nodes)} nodes, {len(synth_graph.edges)} edges.")
        return graph_json_utils.graph_to_json(synth_graph)

    def identify_synth_graph_routes(self, graph_json: SynthGraphJson) -> Tuple[List[SynthRoute], List[SynthRoute], List[CombinationGraph]]:
        if self.enable_profiler:
            logger.debug(
                "Enabling profiler for identify_synth_graph_routes...")
            self.profiler.enable()

        try:
            synth_graph = graph_json_utils.json_to_graph(
                graph_json.model_dump())
        except Exception:
            logger.error("Error parsing synthesis graph JSON:", exc_info=True)
            raise SynthGraphParsingException(
                "Error parsing synthesis graph JSON")

        target_nodes = [
            (node, attrs) for node, attrs in synth_graph.nodes(data=True) if (attrs.get("node_type", "").lower() == "substance" and attrs.get("srole") == "tm")
        ]

        if len(target_nodes) == 0:
            raise ValueError(
                "No target nodes found with node_type 'substance' and srole 'tm'.")
        elif len(target_nodes) > 1:
            raise ValueError(
                f"More than one target node found: {target_nodes}")

        individual_synthesis_routes, individual_synthesis_routes_candidates, combination_graphs = graph_utils.identify_individual_synthesis_routes(
            synth_graph, target_nodes[0][1]["node_id"]
        )
        routes = [SynthRoute(**route) for route in individual_synthesis_routes]
        individual_synthesis_routes_candidates_mapped = [SynthRoute(
            **route) for route in individual_synthesis_routes_candidates]

        combination_graphs_mapped: List[CombinationGraph] = []
        for idx, combination_graph in enumerate(combination_graphs, start=1):
            combination_graphs_mapped.append(CombinationGraph(
                route_index=idx,
                route_candidate=combination_graph
            ))

        if self.enable_profiler:
            logger.debug(
                "Stopping profiler for identify_synth_graph_routes...")
            self.profiler.disable()
            ps = pstats.Stats(self.profiler).sort_stats("cumulative")
            logger.debug(f"Profiling stats:\n{ps}")

        return routes, individual_synthesis_routes_candidates_mapped, combination_graphs_mapped

    def _fetch_synthesis_graph(self, synth_search: SynthGraphSearch) -> SynthGraph:
        """
        Fetches synthesis graph for the given search parameters.

        Args:
            synth_search (SynthGraphSearch): Search parameters to fetch synthesis graph.

        Returns:
            SynthGraph: Synthesis graph for the given search parameters.
        """
        synth_graph_time_start = time.time()

        # Fetch synthesis graph
        synth_graph = self.data_adapter.fetch_synthesis_graph(synth_search)

        # Remove imcompatiable reactions
        synth_graph.synthesis_graph = graph_utils.remove_incompatible_reactions(
            synth_graph.synthesis_graph,
            synth_graph.target_molecule_node_id
        )

        # Assign inventory roles to substances
        self._assign_substance_inventory_roles(
            synth_graph, target_molecule_inchikey=synth_search.target_molecule_inchikey)

        synth_graph.summary = TopologyInfo(
            num_nodes=len(synth_graph.synthesis_graph.nodes()),
            num_edges=len(synth_graph.synthesis_graph.edges()),
        )

        synth_graph_time = time.time() - synth_graph_time_start
        synth_graph.time_info = TimeInfo(synth_graph_time=synth_graph_time)

        # Compute availability info if requested
        if synth_search.include_availability_info:
            availability_info = self._compute_availability_info(synth_graph)
            synth_graph.availability = availability_info

        # Annotate reactions if requested
        if synth_search.annotate_reactions:
            reaction_annotations = self._annotate_reactions(synth_graph)
            synth_graph.synthesis_graph = self._update_graph_with_reaction_annotations(
                synth_graph.synthesis_graph, reaction_annotations)
            synth_graph = self.update_reaction_provenance(synth_graph)
            synth_graph = self.update_reaction_validation(synth_graph)

        return synth_graph

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
        if synthesis_graph_json is None:
            synth_graph = self._fetch_synthesis_graph(route_search)
        else:
            synth_graph = self._parse_synthesis_graph(synthesis_graph_json)

        enumerate_routes_time_start = time.time()

        # Identify synthesis routes
        routes, route_candidates, combination_graphs = self._identify_synth_graph_routes_from_SynthGraph(
            synth_graph, route_search.include_route_candidates, route_search.include_combination_graphs
        )

        enumerate_routes_time = time.time() - enumerate_routes_time_start
        topology_info = TopologyInfo(num_routes=len(routes), num_combination_graphs=len(
            combination_graphs), num_route_candidates=len(route_candidates))

        G_dag = synth_graph.synthesis_graph.copy()

        individual_synthesis_routes: List[SynthRoute] = routes

        aggregate_yield_time_start = time.time()

        # Get the predicted yields for the synthesis graph
        yields = graph_utils.grab_yield_predicted(G_dag)
        projected_routes: List[ProjectedRoute] = []
        # Project the synthesis graph to a reaction network and annotate it with the predicted yields
        for route in individual_synthesis_routes:
            idx = route.route_index
            if idx is None:
                raise ValueError("Route index not found in route")
            # Project the synthesis tree to a reaction network
            projected_route: DiGraph = graph_utils.project_to_reaction_network(
                route.route_candidate)
            # Annotate the projected route with the predicted yields
            graph_utils.annotate_projected_route_with_yield_info(
                projected_route, yields)
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

        final_projected_routes: List[ProjectedRoute] = self._top_n_yield_routes(
            projected_routes, top_n)

        projected_route_indices = {
            pr.route_index for pr in final_projected_routes}

        # Filter out routes in synth_graph.routes that don't correspond to the selected projected routes
        routes = [
            route for route in routes if route.route_index in projected_route_indices]
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
                    route.aggregated_yield = aggregated_yield
                    break  # Once you find the route, break out of the loop

        # Create the final TopNYieldSynthRoutesResult object
        top_n_results: TopNYieldSynthRoutesResult = TopNYieldSynthRoutesResult(
            search_params=route_search,
            synth_graph=synth_graph,
            artifact_free_graph=SynthView(
                view=synth_graph.synthesis_graph.copy()),
            combination_graphs=combination_graphs,
            route_candidates=SynthRouteContainer(subgraphs=route_candidates),
            routes=routes,
            time_info=TimeInfo(
                synth_graph_time=synth_graph.time_info.synth_graph_time if synth_graph.time_info is not None else None,
                enumerate_routes_time=enumerate_routes_time,
                aggregated_yield_time=aggregate_yield_time,
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
        substances = [node_id for node_id, node_data in synth_graph.synthesis_graph.nodes(
            data=True) if node_data.get("node_type") == "substance"]

        for substance in substances:
            inchikey = synth_graph.synthesis_graph.nodes[substance]["inchikey"]
            in_degree = synth_graph.synthesis_graph.in_degree(substance)
            if inchikey == target_molecule_inchikey:
                srole = "tm"
            # elif self.inventory_adapter.in_inventory_by_inchikey(inchikey):
            #     substance.is_in_inventory = True
            #     srole = "sm"
            elif in_degree == 0 and (synth_graph.search_params and synth_graph.search_params.leaves_as_sm):
                srole = "sm"
            else:
                srole = "im"
            synth_graph.synthesis_graph.nodes[substance]["srole"] = srole

    def _identify_synth_graph_routes_from_SynthGraph(
        self, synth_graph: SynthGraph, include_route_candidates: bool, include_combination_graphs: bool
    ) -> Tuple[List[SynthRoute], List[SynthRoute], List[CombinationGraph]]:
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
            route_candidates = [SynthRoute(
                **route) for route in individual_synthesis_routes_candidates]
        if include_combination_graphs:
            combination_maps: List[CombinationGraph] = []
            for idx, combination_graph in enumerate(combination_graphs, start=1):
                combination_maps.append(CombinationGraph(
                    route_index=idx,
                    route_candidate=combination_graph
                ))
        else:
            combination_maps = []
        return routes, route_candidates, combination_maps

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

    def _parse_synthesis_graph(self, synthesis_graph_json: Any) -> SynthGraph:
        """
        Parses synthesis_graph_json into a SynthGraph object with validation steps.

        Raises:
            SynthGraphParsingException: if any step of parsing or validation fails.

        Returns:
            SynthGraph: A validated synthesis graph with metadata.
        """
        # Step 1: Parse JSON to NetworkX DiGraph
        try:
            synth_graph_digraph = graph_json_utils.json_to_graph(
                synthesis_graph_json)
        except Exception as e:
            logger.error(
                "Failed to parse synthesis_graph_json into graph", exc_info=True)
            raise SynthGraphParsingException(
                "Failed to parse synthesis_graph_json into graph") from e

        # Step 2: Validate DAG structure
        try:
            if not graph_utils.is_directed_acyclic_graph(synth_graph_digraph):
                raise SynthGraphParsingException(
                    "Synthesis graph is not a DAG.")
        except Exception as e:
            logger.error("Error validating DAG structure", exc_info=True)
            raise SynthGraphParsingException(
                "Graph structure validation failed") from e

        # Step 3: Identify target node
        try:
            target_nodes = [
                (node, attrs)
                for node, attrs in synth_graph_digraph.nodes(data=True)
                if attrs.get("node_type", "").lower() == "substance" and attrs.get("srole") == "tm"
            ]
            if len(target_nodes) == 0:
                raise ValueError(
                    "No target nodes found with node_type 'substance' and srole 'tm'.")
            elif len(target_nodes) > 1:
                raise ValueError(
                    f"More than one target node found: {target_nodes}")
            target_node_id = target_nodes[0][0]
        except Exception as e:
            logger.error("Error identifying target node", exc_info=True)
            raise SynthGraphParsingException(
                "Failed to identify a valid target node") from e

        # Step 4: Build SynthGraph with metadata
        try:
            parse_time = time.time()
            synth_graph = SynthGraph(
                synthesis_graph=synth_graph_digraph,
                target_molecule_node_id=target_node_id,
                summary=TopologyInfo(
                    num_nodes=len(synth_graph_digraph.nodes()),
                    num_edges=len(synth_graph_digraph.edges()),
                ),
                time_info=TimeInfo(synth_graph_time=time.time() - parse_time)
            )
            return synth_graph
        except Exception as e:
            logger.error("Error constructing SynthGraph object", exc_info=True)
            raise SynthGraphParsingException(
                "Failed to construct SynthGraph") from e

    def aggregate_route_yield(self, synth_route: SynthRoute) -> SynthRoute:
        """
        Aggregates the yield of a given synthesis route
        """
        # Get the predicted yields for the synthesis graph
        yields = graph_utils.grab_yield_predicted(synth_route.route_candidate)
        logger.debug(f"Yields: {yields}")

        # Project route to a reaction network
        projected_route: DiGraph = graph_utils.project_to_reaction_network(
            synth_route.route_candidate)
        # Ensure projected route has yield information
        graph_utils.annotate_projected_route_with_yield_info(
            projected_route, yields)
        # Calculate ag yield
        route_aggregated_yield = ay.aggregate_yields(projected_route)
        # Return route with ag yield
        synth_route.aggregated_yield = route_aggregated_yield
        return synth_route

    def _annotate_reactions(self, synth_graph: SynthGraph) -> Dict[str, Dict[str, str]]:
        """
        Annotates reaction nodes in the synthesis graph with rxclass and rxname.
        Args:
            synth_graph (SynthGraph): Synthesis graph containing reaction nodes.
        """
        reaction_annotations = {}
        for node_id, node_data in synth_graph.synthesis_graph.nodes(data=True):
            if node_data.get("node_type").lower() == "reaction":
                rxid = node_data.get("rxid")
                annotation = self.data_adapter.query_rxid_rxname_rxclass(rxid)
                reaction_annotations[rxid] = annotation

        return reaction_annotations

    def _compute_availability_info(self, synth_graph: SynthGraph) -> List[Availability]:
        """
        Computes availability information for substances in the synthesis graph.
        Args:
            synth_graph (SynthGraph): Synthesis graph containing substance nodes.
        Returns:
            List[Availability]: Availability information for each substance node.
        """
        availability_info: List[Availability] = []
        for node_id, node_data in synth_graph.synthesis_graph.nodes(data=True):
            if node_data.get("node_type") == "substance":
                availability_info.append(
                    # TEMPORARILY SET ALL SM AS AVAILABLE
                    # TODO : Hook into substance availability and ASI inventory
                    Availability(
                        inchikey=node_data.get("inchikey"),
                        inventory=Inventory(available=node_data.get(
                            "srole") == "sm", locations=[]),
                        commercial_availability=CommercialAvailability(
                            available=False, vendors=[]),
                    )
                )
        return availability_info

    def update_reaction_provenance(self, synth_graph: SynthGraph) -> SynthGraph:
        """
        Iterates over the SynthGraph, retrieves source information for each reaction node,
        and updates the patents part of the provenance for reaction nodes.
        Args:
            synth_graph (SynthGraph): The synthesis graph containing reaction nodes.
        Returns:
            SynthGraph: The updated synthesis graph with provenance information.
        """
        for node_id, node_data in synth_graph.synthesis_graph.nodes(data=True):
            if node_data["node_type"].lower() == "reaction":
                rxid = node_data.get("rxid")
                source_info = self.get_reaction_source_info(rxid)
                node_data["original_rxsmiles"] = source_info["rxsmiles_original"]
                # because only one patent source for now
                node_data["patents"] = [source_info["patents"]]
                # because only one patent source for now
                node_data["patent_paragraph_nums"] = [
                    source_info["paragraphs"]]
        return synth_graph

    def update_reaction_validation(self, synth_graph: SynthGraph) -> SynthGraph:
        """
        Iterates over the SynthGraph, checks the rxname of each reaction node,
        and sets is_rxname_recognized to True if rxname is not 'Unknown' and not an empty string.
        Args:
            synth_graph (SynthGraph): The synthesis graph containing reaction nodes.
        Returns:
            SynthGraph: The updated synthesis graph with validation information.
        """
        for node_id, node_data in synth_graph.synthesis_graph.nodes(data=True):
            if node_data["node_type"].lower() == "reaction":
                rxname = node_data.get("rxname", "")
                if rxname and rxname not in ["Unknown", "Unrecognized"]:
                    node_data["is_rxname_recognized"] = True
        return synth_graph

    def get_reaction_source_info(self, rxid: str) -> Dict[str, Any]:
        """
        Retrieves source information for a reaction, including patents, rxsmiles_original, and other metadata.
        Args:
            rxid (str): The reaction ID.
        Returns:
            Dict[str, Any]: A dictionary containing source information for the reaction.
        """
        reaction = self.data_adapter.get_reaction_by_rxid(rxid)
        patent_number = reaction["PatentNumber"]
        paragraph_number = reaction["ParagraphNum"]
        rxsmiles_original = reaction.get("rxsmiles_original", None)

        return {"patents": patent_number, "paragraphs": paragraph_number, "rxsmiles_original": rxsmiles_original}

    def _update_graph_with_reaction_annotations(self, graph: DiGraph, reaction_annotations: Dict[str, Dict[str, str]]) -> DiGraph:
        """
        Updates the rxclass and rxname fields for reaction nodes in the given graph.
        Args:
            graph (DiGraph): The graph to update.
            reaction_annotations (Dict[str, Dict[str, str]]): A dictionary mapping node IDs to rxclass and rxname.
        """
        for rxid, annotations in reaction_annotations.items():
            for node_id, node_data in graph.nodes(data=True):
                if node_data.get("node_label") == rxid:
                    graph.nodes[node_id]["rxclass"] = annotations["rxclass"]
                    graph.nodes[node_id]["rxname"] = annotations["rxname"]
                    break

        return graph
