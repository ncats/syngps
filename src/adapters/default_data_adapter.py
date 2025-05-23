import json
from typing import Any, Dict, List, Optional, Set

import pandas as pd
from networkx import DiGraph, all_simple_paths
from adapters import AbstractDataAdapter
from errors import (
    CsvLoadError,
    GraphLoadError,
    InchikeyNotFoundError,
    MultipleInchikeyFoundError,
    ReactionNotFoundError,
)
from models import SynthGraph, SynthGraphSearch
import logging

logger = logging.getLogger(__name__)

class DefaultDataAdapter(AbstractDataAdapter):
    """
    Default Data Adapter Class. This class implements the AbstractDataAdapter interface.
    It is used to interact with the data sources in a csv and jsonl format.
    """

    def __init__(
        self,
        csv_file_path: str = "./reactions_data.csv",
        jsonl_file_path: str = "./xs_small_graph.jsonl",
    ):
        """
        Constructor for the DefaultDataAdapter class.

        Args:
            csv_file_path (str): Path to the CSV file containing reaction data.
            jsonl_file_path (str): Path to the JSONL file containing graph data.

        Raises:
            CsvLoadError: If the CSV file fails to load.
            GraphLoadError: If the JSONL file fails to load.
        """
        # Load CSV data
        try:
            logger.debug(f"Loading reaction CSV file from {csv_file_path}")
            with open(csv_file_path, "r") as file:
                self.data = pd.read_csv(file, dtype={"rxid": str})
        except Exception as e:
            logger.error(f"Failed to load CSV file at {csv_file_path}", exc_info=True)
            raise CsvLoadError(f"Failed to load CSV file: {e}")

        # Load JSONL data into a NetworkX graph
        self.graph = DiGraph()
        try:
            logger.debug(f"Loading graph JSON file from {jsonl_file_path}")
            with open(jsonl_file_path, "r") as file:
                for line in file:
                    data = json.loads(line)

                    if data["type"] == "node":
                        # Add node with properties
                        node_id = data["id"]
                        properties = data.get("properties", {})
                        if properties.get("node_type", None) is None:
                            labels = data.get("labels", ["unknown"])
                            properties["node_type"] = labels[0].lower()
                        self.graph.add_node(node_id, **properties)

                    elif data["type"] == "relationship":
                        # Add edge with properties
                        start_node = data["start"]["id"]
                        end_node = data["end"]["id"]
                        properties = data.get("properties", {})
                        edge_label: str = data["label"]
                        properties["edge_type"] = edge_label.lower()
                        # Including the edge label as a property
                        self.graph.add_edge(start_node, end_node, label=edge_label, **properties)
        except Exception as e:
            logger.error(f"Failed to load graph JSONL file at {jsonl_file_path}", exc_info=True)
            raise GraphLoadError(f"Failed to load Graph JSONL file: {e}")

    # TODO : Create Reaction model
    def get_reaction_by_rxid(self, rxid: str) -> Dict[str, Any]:
        """
        Method to get reaction for a given reaction ID.

        Args:
            rxid (str): Reaction ID.

        Returns:
            Optional[Dict[str, Any]]: Reaction details if found, None otherwise.

        Raises:
            ReactionNotFoundError: If the reaction ID is not found.
        """
        result = self.data[self.data["rxid"] == rxid]
        if not result.empty:
            return result.iloc[0]
        else:
            logger.debug(f"Reaction with rxid {rxid} not found, raising ReactionNotFoundError")
            raise ReactionNotFoundError(rxid)

    def verify_connection(self) -> bool:
        """
        Method to verify connection to the knowledge base.

        Returns:
            bool: True if connection is successful, False otherwise.
        """
        return True

    def close(self):
        """
        No action required to close the connection to the knowledge base.
        """
        pass

    def get_rxsmiles_extended(self, rxid: str) -> str:
        """
        Method to get extended rxsmiles for a given reaction ID.

        Args:
            rxid (str): Reaction ID.

        Returns:
            str: Extended rxsmiles if found, None otherwise.

        Raises:
            ReactionNotFoundError: If the reaction ID is not found.
        """
        result = self.data[self.data["rxid"] == rxid]["rxsmiles"]
        if not result.empty:
            return result.iloc[0]
        else:
            logger.debug(f"Reaction with rxid {rxid} not found, raising ReactionNotFoundError")
            raise ReactionNotFoundError(rxid)

    def get_multiple_rxsmiles_extended(self, rxids: List[str]) -> Dict[str, Optional[str]]:
        """
        Method to get multiple extended rxsmiles for a list of reaction IDs. Value is None if rxid is not found.

        Args:
            rxids (List[str]): List of reaction IDs.

        Returns:
            Dict[str, Optional[str]]: Dictionary of rxid to extended rxsmiles.
        """
        rxsmiles_extended_dict: Dict[str, Optional[str]] = {}
        for rxid in rxids:
            result = self.data[self.data["rxid"] == rxid]["rxsmiles"]
            rxsmiles_extended_dict[rxid] = result.iloc[0] if not result.empty else None
        return rxsmiles_extended_dict

    def find_target_molecule_node(self, target_molecule: str) -> str:
        """
        Method to find the target molecule node.

        Args:
            target_molecule (str): Target molecule.

        Returns:
            str: Node ID of the target molecule.

        Raises:
            InchikeyNotFoundError: If the target molecule inchikey is not found.
            MultipleInchikeyFoundError: If multiple nodes with the same inchikey are found.
        """
        target_nodes = [n for n, attr in self.graph.nodes(data=True) if attr.get("inchikey") == target_molecule and attr.get("node_type") == "substance"]
        if len(target_nodes) == 0:
            logger.debug(f"Target molecule with inchikey {target_molecule} not found, raising InchikeyNotFoundError")
            raise InchikeyNotFoundError(target_molecule)
        elif len(target_nodes) > 1:
            logger.debug(f"Multiple nodes found for target molecule with inchikey {target_molecule}, raising MultipleInchikeyFoundError")
            raise MultipleInchikeyFoundError(target_molecule)
        else:
            return target_nodes[0]

    def fetch_ambiguous_reactions(self) -> Set[str]:
        """
        Method to find ambiguous reactions in the knowledge base.

        Returns:
            Set[str]: A set of distinct reaction IDs (r.rxid)
        """
        reaction_ids: Set[str] = set()

        for node_id, node_properties in self.graph.nodes(data=True):
            if node_properties.get("node_type") == "reaction":
                products: Set[int] = set()
                inputs: Set[int] = set()

                # Loop over incoming edges
                for start_node_id, _, edge_data in self.graph.in_edges(node_id, data=True):
                    if edge_data.get("edge_type").upper() == "REACTANT_OF" or edge_data.get("edge_type").upper() == "REAGENT_OF":
                        inputs.add(start_node_id)

                # Loop over outgoing edges
                for _, end_node_id, edge_data in self.graph.out_edges(node_id, data=True):
                    if edge_data.get("edge_type").upper() == "PRODUCT_OF":
                        products.add(end_node_id)

                # Check for overlap between products and inputs
                if len(products.intersection(inputs)) > 0:
                    reaction_ids.add(node_properties.get("rxid"))

        return reaction_ids

    def fetch_multiproduct_reactions(self) -> Set[str]:
        """
        Method to find multiproduct reactions in the knowledge base.

        Returns:
            Set[str]: Set of multiproduct reactions.
        """
        reaction_ids = set()

        for node in self.graph.nodes(data=True):
            if node[1].get("node_type") == "reaction":
                product_count = 0

                # Count 'PRODUCT_OF' relationships for this reaction
                for _, _, edge_data in self.graph.edges(node[0], data=True):
                    if edge_data.get("edge_type").upper() == "PRODUCT_OF".upper():
                        product_count += 1

                # If the reaction has more than one product, add its rxid to the set
                if product_count > 1:
                    reaction_ids.add(node[1].get("rxid"))

        return reaction_ids

    def fetch_synthesis_graph(self, search_params: SynthGraphSearch) -> SynthGraph:
        """
        Method to fetch the synthesis graph based on search parameters.

        Args:
            search_params (SynthGraphSearch): Search parameters.

        Returns:
            SynthGraph: Synthesis graph.

        Raises:
            InchikeyNotFoundError: If the target molecule inchikey is not found.
        """
        logger.debug(f"Fetching synthesis graph with search params: {search_params}")
        target_node = self.find_target_molecule_node(search_params.target_molecule_inchikey)
        G = self._nx_find_reactions(
            self.graph.copy(),
            search_params.target_molecule_inchikey,
            search_params.search_depth,
            search_params.query_type,
        )
        logger.debug(f"Synthesis graph fetched: {G}")

        for node_id in G.nodes:
            # TODO : Check if this is needed for networkx graphs
            # Set node_id for each node in the graph (if necessary)
            G.nodes[node_id]["node_id"] = node_id

        if target_node is None:
            logger.debug(f"Target molecule with inchikey {search_params.target_molecule_inchikey} not found, raising InchikeyNotFoundError")
            raise InchikeyNotFoundError(search_params.target_molecule_inchikey)

        synth_graph = SynthGraph(
            target_molecule_node_id=target_node,
            synthesis_graph=G.copy(),
            search_params=search_params,
        )

        return synth_graph

    def get_data_counts(self):
        raise NotImplementedError("Method not implemented for default data adapter yet.")

    def get_substance_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        return {"inchikey": inchikey, "nsinchikey": inchikey, "canonical_smiles": "NA"}

    # ====== HELPERS ======
    def _nx_find_reactions(
        self,
        graph,
        target_molecule,
        search_depth,
        query_type="shortest_path",
    ):
        """
        Method to find reactions in a NetworkX graph.

        Args:
            graph (DiGraph): NetworkX graph.
            target_molecule (str): Target molecule inchikey.
            search_depth (int): Search depth.
            query_type (str): Query type.

        Returns:
            DiGraph: Subgraph containing the reactions.
        """
        # Sets to keep track of nodes and edges in valid paths
        subgraph_nodes = set()
        subgraph_edges = set()
        if query_type == "shortest_path":
            target_nodes = [n for n, attr in graph.nodes(data=True) if attr.get("inchikey") == target_molecule and attr.get("node_type") == "substance"]

            if len(target_nodes) == 0:
                raise InchikeyNotFoundError(target_molecule)
            elif len(target_nodes) > 1:
                raise MultipleInchikeyFoundError(target_molecule)

            for target in target_nodes:
                for node in graph.nodes:
                    if graph.nodes[node].get("node_type") == "reaction":
                        for path in all_simple_paths(graph, source=node, target=target, cutoff=search_depth):
                            # Add nodes and edges from the valid path to the sets
                            subgraph_nodes.update(path)
                            # Assert all incoming and outgoing edges for each reaction node are in the sets
                            for n in path:
                                if graph.nodes[n].get("node_type") == "reaction":
                                    for u, v in graph.in_edges(n):
                                        subgraph_edges.add((u, v))
                                        subgraph_nodes.add(u)
                                    for u, v in graph.out_edges(n):
                                        subgraph_edges.add((u, v))
                                        subgraph_nodes.add(v)

            subgraph = DiGraph()
            # Create subgraph from the collected nodes and edges
            subgraph.add_nodes_from((n, graph.nodes[n]) for n in subgraph_nodes)
            subgraph.add_edges_from((u, v, graph[u][v]) for u, v in subgraph_edges if u in subgraph_nodes and v in subgraph_nodes)
        elif query_type == "full_graph":
            raise NotImplementedError("'full_graph' query type not implemented for default data adapter yet.")
        else:
            raise Exception("[ERROR] Invalid query type. Source: networkx_graph_utils.py, fn nx_find_reactions().")

        return subgraph
