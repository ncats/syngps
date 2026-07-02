import os
from typing import Any, Dict, List, Optional, Set
import pickle
from networkx import DiGraph, all_simple_paths, NetworkXNoPath, node_link_graph
from syngps.adapters import AbstractDataAdapter
from syngps.errors import (
    GraphLoadError,
    InchikeyNotFoundError,
    MultipleInchikeyFoundError,
    MultipleReactionFoundError,
    ReactionNotFoundError,
)
from syngps.models import SynthGraph, SynthGraphSearch
from syngps.utils import safe_node_link_data
import logging

logger = logging.getLogger(__name__)

DEFAULT_PKL = os.path.join(os.path.dirname(__file__), "aicp_1k_rxn_subset.pkl")


class NetworkxDataAdapter(AbstractDataAdapter):
    """
    Networkx Data Adapter. Loads a networkx graph from a pickled file or JSON object.
    """

    networkx_graph: DiGraph

    def __init__(
        self,
        aicp_pkl_path: Optional[str] = DEFAULT_PKL,
        json_data: Optional[dict] = None
    ):
        if json_data:
            try:
                logger.debug("Loading networkx graph from JSON input.")
                graph = node_link_graph(json_data, directed=True)

                if not isinstance(graph, DiGraph):
                    raise GraphLoadError(
                        f"Loaded object is not a DiGraph: {type(graph).__name__}"
                    )

                self.networkx_graph = graph
                logger.debug("Successfully loaded graph from JSON.")

            except Exception as e:
                logger.error(
                    "Failed to load graph from JSON data.", exc_info=True)
                raise GraphLoadError(
                    f"Failed to load networkx graph from JSON: {e}")
        else:
            if not aicp_pkl_path or not os.path.exists(aicp_pkl_path):
                raise FileNotFoundError(
                    f"Pickled data file not found: {aicp_pkl_path}")
            self.aicp_pkl_path = aicp_pkl_path

            try:
                with open(aicp_pkl_path, "rb") as file:
                    logger.debug(
                        f"Loading networkx graph from {aicp_pkl_path}")
                    graph = pickle.load(file)

                    if not isinstance(graph, DiGraph):
                        raise GraphLoadError(
                            f"Loaded object is not a NetworkX graph: {type(graph).__name__}"
                        )

                    self.networkx_graph = graph
                    logger.debug("Successfully loaded graph from PKL.")
            except GraphLoadError:
                raise
            except Exception as e:
                logger.error(
                    f"Failed to load networkx graph from {aicp_pkl_path}", exc_info=True)
                raise GraphLoadError(f"Failed to load networkx graph: {e}")

    def export_networkx_graph_to_json(self):
        """
        Export the networkx graph to a JSON object.
        """
        return safe_node_link_data(self.networkx_graph)

    def verify_connection(self):
        """
        Verify connection to networkx graph, ensuring it exists.
        """
        return bool(self.networkx_graph)

    def close(self):
        """
        Delete the networkx graph
        """
        pass

    def get_data_counts(self):
        """
        Get data counts for the networkx graph.
        """
        return {
            "total_nodes": len(self.networkx_graph.nodes),
            "total_edges": len(self.networkx_graph.edges),
            "total_reactions": self._count_node_types("reaction"),
            "total_substances": self._count_node_types("substance"),
            "total_product_of": self._count_edge_types("PRODUCT_OF"),
            "total_reactant_of": self._count_edge_types("REACTANT_OF"),
            "total_reagent_of": self._count_edge_types("REAGENT_OF"),
        }

    def get_reaction_by_rxid(self, rxid: str) -> Dict[str, Any]:
        """
        Get a reaction node by its rxid.

        Args:
            rxid (str): Reaction ID to search for.

        Returns:
            Dict[str, Any]: Node attributes of the matching reaction.

        Raises:
            ReactionNotFoundError: If no matching reaction is found.
            MultipleInchikeyFoundError: If more than one matching node is found.
        """
        matches = [
            (node_id, data)
            for node_id, data in self.networkx_graph.nodes(data=True)
            if data.get("node_type").lower() == "reaction" and data.get("rxid") == rxid
        ]

        if not matches:
            raise ReactionNotFoundError(
                f"No reaction found with rxid '{rxid}'.")

        if len(matches) > 1:
            raise MultipleReactionFoundError(
                f"Multiple reactions found with rxid '{rxid}'.")

        node_id, node_data = matches[0]
        # Include node_id in the returned dict for reference
        return {"node_id": node_id, **node_data}

    def get_rxsmiles_extended(self, rxid: str) -> str:
        """
        Method to get extended rxsmiles for a given reaction ID.

        Args:
            rxid (str): Reaction ID.

        Returns:
            str: Extended rxsmiles if found, None otherwise.
        """
        return self.get_reaction_by_rxid(rxid)["rxsmiles"]

    def get_multiple_rxsmiles_extended(self, rxids: List[str]) -> Dict[str, Optional[str]]:
        """
        Method to get multiple extended rxsmiles for a list of reaction IDs.

        Args:
            rxids (List[str]): List of reaction IDs.

        Returns:
            Dict[str, Optional[str]]: Dictionary of rxid to extended rxsmiles.
        """
        return {rxid: self.get_rxsmiles_extended(rxid) for rxid in rxids}

    def get_substance_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Method to get substance for a given InChIKey.

        Args:
            inchikey (str): InChIKey.

        Returns:
            Optional[Dict[str, Any]]: Substance details if found, None otherwise.
        """
        matches = [
            (node_id, data)
            for node_id, data in self.networkx_graph.nodes(data=True)
            if data.get("node_type").lower() == "substance" and data.get("inchikey") == inchikey
        ]

        if not matches:
            raise InchikeyNotFoundError(
                f"No node found with inchikey '{inchikey}'.")

        if len(matches) > 1:
            raise MultipleInchikeyFoundError(
                f"Multiple nodes found with inchikey '{inchikey}'.")

        return matches[0][1]

    def find_target_molecule_node(self, target_molecule: str) -> str:
        """
        Method to find the target molecule node.

        Args:
            target_molecule (str): Inchikey of the target molecule.

        Returns:
            str: Node ID of the target molecule.
        """
        matches = [
            (node_id, data)
            for node_id, data in self.networkx_graph.nodes(data=True)
            if data.get("node_type").lower() == "substance" and data.get("inchikey") == target_molecule
        ]

        if not matches:
            raise InchikeyNotFoundError(
                f"No node found with inchikey '{target_molecule}'.")

        if len(matches) > 1:
            raise MultipleInchikeyFoundError(
                f"Multiple nodes found with inchikey '{target_molecule}'.")

        return matches[0][0]

    def query_rxid_rxname_rxclass(self, rxid: str) -> Dict[str, str]:
        """
        Queries the NetworkX graph for a reaction node and its associated RXName and RXClass nodes.

        Args:
            rxid (str): The reaction ID to query.

        Returns:
            Dict[str, str]: A dictionary containing 'rxname' and 'rxclass'.
        """
        # Step 1: Find the Reaction node with the matching rxid
        reaction_node = None
        for node_id, data in self.networkx_graph.nodes(data=True):
            if data.get("node_type").lower() == "reaction" and data.get("rxid") == rxid:
                reaction_node = node_id
                break

        if not reaction_node:
            return {"rxname": "", "rxclass": ""}

        # Step 2: Find the RXName node connected to the Reaction via RXNM2R
        rxname_node = None
        for source, target, edge_data in self.networkx_graph.in_edges(reaction_node, data=True):
            if (
                edge_data.get("edge_type") == "RXNM2R" and
                self.networkx_graph.nodes[source].get("node_type") == "RXName"
            ):
                rxname_node = source
                break

        if not rxname_node:
            return {"rxname": "", "rxclass": ""}

        # Step 3: Find the RXClass node connected to the RXName via RXCL2RXNM
        rxclass_node = None
        for source, target, edge_data in self.networkx_graph.in_edges(rxname_node, data=True):
            if (
                edge_data.get("edge_type") == "RXCL2RXNM" and
                self.networkx_graph.nodes[source].get("node_type") == "RXClass"
            ):
                rxclass_node = source
                break

        # Step 4: Extract values safely
        rxname = self.networkx_graph.nodes[rxname_node].get(
            "rxname", "") if rxname_node else ""
        rxclass = self.networkx_graph.nodes[rxclass_node].get(
            "rxclass", "") if rxclass_node else ""

        return {
            "rxname": rxname,
            "rxclass": rxclass,
        }

    def fetch_ambiguous_reactions(self) -> Set[str]:
        """
        Method to find ambiguous reactions where a product is also a reactant or reagent
        in the same reaction.

        Returns:
            Set[str]: Set of ambiguous reaction rxids.
        """
        ambiguous_rxids = set()

        for node_id, node_data in self.networkx_graph.nodes(data=True):
            # Focus only on Reaction nodes
            if node_data.get("node_type").lower() != "reaction":
                continue

            rxid = node_data.get("rxid")
            if not rxid:
                continue

            # Track substances by edge role
            product_substances = set()
            input_substances = set()

            for source, target, edge_data in self.networkx_graph.in_edges(node_id, data=True):
                if edge_data.get("edge_type") in ("REACTANT_OF", "REAGENT_OF"):
                    input_substances.add(source)

            for source, target, edge_data in self.networkx_graph.out_edges(node_id, data=True):
                if edge_data.get("edge_type") == "PRODUCT_OF":
                    product_substances.add(target)

            # Check for overlap
            if product_substances & input_substances:
                ambiguous_rxids.add(rxid)

        return ambiguous_rxids

    def fetch_multiproduct_reactions(self) -> Set[str]:
        """
        Method to find multiproduct reactions in the knowledge base.

        Returns:
            Set[str]: Set of multiproduct reaction rxids.
        """
        multiproduct_rxids = set()

        for node_id, node_data in self.networkx_graph.nodes(data=True):
            # Only consider Reaction nodes
            if node_data.get("node_type").lower() != "reaction":
                continue

            rxid = node_data.get("rxid")
            if not rxid:
                continue

            # Find all PRODUCT_OF edges from this reaction
            product_targets = {
                target for _, target, edge_data in self.networkx_graph.out_edges(node_id, data=True)
                if edge_data.get("edge_type") == "PRODUCT_OF"
            }

            if len(product_targets) > 1:
                multiproduct_rxids.add(rxid)

        return multiproduct_rxids

    def fetch_synthesis_graph(self, search_params):
        """
        Method to fetch the synthesis graph based on search parameters.

        Args:
            search_params (SynthGraphSearch): Search parameters.

        Returns:
            SynthGraph: Synthesis graph.
        """
        return self._build_synthesis_graph_from_networkx(self.networkx_graph, search_params)

    def get_reaction_products(self, rxid: str) -> List[Dict[str, Any]]:
        """
        Retrieves the products of a reaction node from the NetworkX graph by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve the products for.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the product nodes for the reaction.
        """
        products: List[Dict[str, Any]] = []

        # Find the reaction node
        reaction_node = next(
            (n for n, d in self.networkx_graph.nodes(data=True)
             if d.get("node_type") == "reaction" and d.get("rxid") == rxid),
            None
        )
        if not reaction_node:
            return products

        # Find all PRODUCT_OF edges from reaction to product substances
        for _, target, edge_data in self.networkx_graph.out_edges(reaction_node, data=True):
            if edge_data.get("edge_type") == "PRODUCT_OF":
                target_data = self.networkx_graph.nodes[target]
                inchikey = target_data.get("inchikey")
                if inchikey:
                    substance = self.get_substance_by_inchikey(inchikey)
                    if substance:
                        products.append(substance)
                    else:
                        raise InchikeyNotFoundError(
                            f"Substance with inchikey '{inchikey}' not found.")

        return products

    def get_reaction_reactants(self, rxid: str) -> List[Dict[str, Any]]:
        """
        Retrieves the reactants of a reaction node from the NetworkX graph by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve the reactants for.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the reactant nodes for the reaction.
        """
        reactants: List[Dict[str, Any]] = []

        reaction_node = next(
            (n for n, d in self.networkx_graph.nodes(data=True)
             if d.get("node_type") == "reaction" and d.get("rxid") == rxid),
            None
        )
        if not reaction_node:
            return reactants

        for source, _, edge_data in self.networkx_graph.in_edges(reaction_node, data=True):
            if edge_data.get("edge_type") == "REACTANT_OF":
                source_data = self.networkx_graph.nodes[source]
                inchikey = source_data.get("inchikey")
                if inchikey:
                    substance = self.get_substance_by_inchikey(inchikey)
                    if substance:
                        reactants.append(substance)
                    else:
                        raise InchikeyNotFoundError(
                            f"Substance with inchikey '{inchikey}' not found.")

        return reactants

    def get_reaction_reagents(self, rxid: str) -> List[Dict[str, Any]]:
        """
        Retrieves the reagents of a reaction node from the NetworkX graph by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve the reagents for.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the reagent nodes for the reaction.
        """
        reagents: List[Dict[str, Any]] = []

        reaction_node = next(
            (n for n, d in self.networkx_graph.nodes(data=True)
             if d.get("node_type") == "reaction" and d.get("rxid") == rxid),
            None
        )
        if not reaction_node:
            return reagents

        for source, _, edge_data in self.networkx_graph.in_edges(reaction_node, data=True):
            if edge_data.get("edge_type") == "REAGENT_OF":
                source_data = self.networkx_graph.nodes[source]
                inchikey = source_data.get("inchikey")
                if inchikey:
                    substance = self.get_substance_by_inchikey(inchikey)
                    if substance:
                        reagents.append(substance)
                    else:
                        raise InchikeyNotFoundError(
                            f"Substance with inchikey '{inchikey}' not found.")

        return reagents

    ###################
    # HELPER FUNCTIONS
    ###################

    def _count_node_types(self, node_type):
        return len([node for node in self.networkx_graph.nodes if self.networkx_graph.nodes[node]["node_type"].lower() == node_type.lower()])

    def _count_edge_types(self, edge_type):
        return len([edge for edge in self.networkx_graph.edges if self.networkx_graph.edges[edge]["edge_type"].lower() == edge_type.lower()])

    def _build_synthesis_graph_from_networkx(self, networkx_graph: DiGraph, search_params: SynthGraphSearch) -> SynthGraph:
        """
        Builds a synthesis graph from the NetworkX graph based on the given search parameters.
        Args:
            networkx_graph (DiGraph): The NetworkX graph to use for synthesis graph operations.
            search_params (SynthGraphSearch): The search parameters to use for building the synthesis graph.
        Returns:
            SynthGraph: A SynthGraph object containing the synthesis graph and other relevant information.
        """
        logger.info("Building synthesis graph from NetworkX graph.")
        target_inchikey = search_params.target_molecule_inchikey
        search_depth = search_params.search_depth
        valid_rels = {"PRODUCT_OF", "REACTANT_OF", "REAGENT_OF"}

        # Step 1: Find target node by inchikey
        target_node = None
        for node, data in networkx_graph.nodes(data=True):
            if data.get("node_type").lower() == "substance" and data.get("inchikey") == target_inchikey:
                target_node = node
                break
        if not target_node:
            raise ValueError(
                f"Target molecule with inchikey '{target_inchikey}' not found in graph.")

        # Step 2: Find all upstream substance nodes (excluding the target)
        upstream_substances = set()

        def dfs_upstream(node, depth):
            if depth == 0:
                return
            for pred in networkx_graph.predecessors(node):
                edge_data = networkx_graph.get_edge_data(pred, node)
                if edge_data and edge_data.get("edge_type").upper() in {"PRODUCT_OF", "REACTANT_OF"}:
                    pred_data = networkx_graph.nodes[pred]
                    if pred_data.get("node_type").lower() == "substance" and pred_data.get("inchikey") != target_inchikey:
                        upstream_substances.add(pred)
                    dfs_upstream(pred, depth - 1)

        dfs_upstream(target_node, search_depth)

        # Step 3: Collect all valid paths from each upstream substance to the target
        path_nodes = set()
        path_edges = set()

        for source in upstream_substances:
            try:
                for path in all_simple_paths(
                    networkx_graph,
                    source=source,
                    target=target_node,
                    cutoff=search_depth
                ):
                    for i in range(len(path) - 1):
                        u, v = path[i], path[i + 1]
                        edge_data = networkx_graph.get_edge_data(u, v)
                        if edge_data and edge_data.get("edge_type").upper() in {"PRODUCT_OF", "REACTANT_OF"}:
                            path_nodes.update([u, v])
                            path_edges.add((u, v))
            except NetworkXNoPath:
                continue

        # Step 4: Collect additional edges from Reaction nodes (REACTANT_OF, PRODUCT_OF, REAGENT_OF)
        full_nodes = set(path_nodes)
        full_edges = set(path_edges)

        for node in path_nodes:
            data = networkx_graph.nodes[node]
            if data.get("node_type").lower() == "reaction":
                for nbr in networkx_graph.successors(node):
                    edge_data = networkx_graph.get_edge_data(node, nbr)
                    if edge_data and edge_data.get("edge_type").upper() in valid_rels:
                        full_nodes.add(nbr)
                        full_edges.add((node, nbr))
                for nbr in networkx_graph.predecessors(node):
                    edge_data = networkx_graph.get_edge_data(nbr, node)
                    if edge_data and edge_data.get("edge_type").upper() in valid_rels:
                        full_nodes.add(nbr)
                        full_edges.add((nbr, node))

        # Step 5: Construct subgraph
        synthesis_graph = networkx_graph.subgraph(full_nodes).copy()

        # Add only the filtered edges
        synthesis_graph = DiGraph()
        for node in full_nodes:
            synthesis_graph.add_node(node, **networkx_graph.nodes[node])
        for u, v in full_edges:
            synthesis_graph.add_edge(
                u, v, **networkx_graph.get_edge_data(u, v))

        # Step 6: Enrich node and edge attributes
        for node in synthesis_graph.nodes:
            updated_properties = self._extract_node_properties_from_dict(
                synthesis_graph.nodes[node])
            synthesis_graph.nodes[node].update(updated_properties)

        for u, v, edge_data in synthesis_graph.edges(data=True):
            edge_data["start_node"] = u
            edge_data["end_node"] = v
            updated_edge_properties = self._extract_edge_properties_from_dict(
                edge_data)
            synthesis_graph[u][v].update(updated_edge_properties)

        return SynthGraph(
            target_molecule_node_id=target_node,
            synthesis_graph=synthesis_graph,
            search_params=search_params,
        )

    def _extract_node_properties_from_dict(self, node_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extracts the properties of a node from a dictionary.
        Args:
            node_data (Dict[str, Any]): The dictionary containing node properties.
        Returns:
            Dict[str, Any]: A dictionary containing the properties of the node.
        Raises:
            Exception: If the node type is not recognized.
        """
        node_properties = node_data.copy()

        if node_properties.get("node_type", None) is None:
            label = node_properties.get("type", "unknown")
            node_properties["node_type"] = label.lower()

        node_type = node_properties.get("node_type", None).lower()
        if node_type == "substance":
            node_properties["node_label"] = node_properties["inchikey"]
        elif node_type == "reaction":
            node_properties["node_label"] = node_properties["rxid"]
            node_properties["yield_predicted"] = node_properties.get(
                "yield_predicted", 0.0)
            node_properties["yield_score"] = node_properties.get(
                "yield_score", 0.0)
        else:
            logger.debug(f"Node type not recognized: {node_type}")
            raise Exception("[ERROR] Invalid node type.")
        return node_properties

    def _extract_edge_properties_from_dict(self, edge_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extracts the properties of an edge from a dictionary.
        Args:
            edge_data (Dict[str, Any]): The dictionary containing edge properties.
        Returns:
            Dict[str, Any]: A dictionary containing the properties of the edge.
        Raises:
            Exception: If the edge type is not recognized.
        """
        edge_properties = edge_data.copy()
        edge_properties["edge_type"] = edge_properties["edge_type"].lower()
        forward_s2r_edge_types = ["reactant_of", "reagent_of"]
        reverse_s2r_edge_types = ["product_of"]

        if edge_properties.get("edge_type") in forward_s2r_edge_types:
            edge_properties["inchikey"] = edge_properties.get("start_node")
            edge_properties["rxid"] = edge_properties.get("end_node")
            edge_properties["start_node"] = edge_properties.get("start_node")
            edge_properties["end_node"] = edge_properties.get("end_node")
        elif edge_properties.get("edge_type") in reverse_s2r_edge_types:
            edge_properties["inchikey"] = edge_properties.get("end_node")
            edge_properties["rxid"] = edge_properties.get("start_node")
            edge_properties["start_node"] = edge_properties.get("start_node")
            edge_properties["end_node"] = edge_properties.get("end_node")
        else:
            raise Exception("[ERROR] Invalid edge type.")

        return edge_properties
