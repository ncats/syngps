from functools import lru_cache
from typing import Any, Dict, List, Optional, Set

from networkx import DiGraph
from pymongo import MongoClient
from syngps import AbstractDataAdapter
from syngps.errors import (
    InchikeyNotFoundError,
    MultipleInchikeyFoundError,
    ReactionNotFoundError,
    SubstanceNotFoundInSynthGraphError,
)
from syngps.models import SynthGraph, SynthGraphSearch
from app.logging_config import logger
from app.utils import decrypt_value_with_prepended_iv_aes_gcm

# Neo4j driver is used to connect to the graph database
from neo4j import GraphDatabase


class KBAdapter(AbstractDataAdapter):
    """
    A data adapter for interacting with the default AICP knowledge base consisting of Memgraph, MongoDB, and Neo4j databases.
    Extends the AbstractDataAdapter class.

    Attributes:
        memgraph_driver (GraphDatabase.memgraph_driver): The driver for connecting to the Memgraph database.
        mongo_client (MongoClient): The client for connecting to the MongoDB database.
        neo4j_driver (GraphDatabase.driver): The driver for connecting to the Neo4j database.
    """

    def __init__(
        self,
        mongo_uri: str,
        mongo_db_name: str,
        memgraph_uri: str,
        graph_backend: str = "memgraph",
        mongo_user: Optional[str] = None,
        mongo_password: Optional[str] = None,
        memgraph_user: Optional[str] = None,
        memgraph_password: Optional[str] = None,
        neo4j_uri: Optional[str] = None,
        neo4j_user: Optional[str] = None,
        neo4j_password: Optional[str] = None,
        memgraph_conn_encrypted: bool = False,
        neo4j_conn_encrypted: bool = False,
        smiles_encryption_key: Optional[str] = None,
    ):
        """
        Initializes the KBAdapter with connections to mongo, memgraph, and neo4j based on the graph_backend choice.

        Args:
            mongo_uri (str): The URI for the MongoDB database.
            mongo_db_name (str): The name of the MongoDB database.
            memgraph_uri (str): The URI for the Memgraph database.
            neo4j_uri (str): The URI for the Neo4j database.
            graph_backend (str): The graph database backend to use. Can be 'memgraph' or 'neo4j'. Default is 'memgraph'.
            mongo_user (str): The username for the MongoDB database.
            mongo_password (str): The password for the MongoDB database.
            memgraph_user (str): The username for the Memgraph database.
            memgraph_password (str): The password for the Memgraph database.
            neo4j_user (str): The username for the Neo4j database.
            neo4j_password (str): The password for the Neo4j database.
        """
        self.REACTION_COLLECTION_NAME = "reactionCollection"
        self.SUBSTANCE_COLLECTION_NAME = "substanceCollection"
        self.SMILES_ENCRYPTION_KEY = smiles_encryption_key
        self.graph_backend = graph_backend

        # Set up connection to MongoDB database
        if mongo_user and mongo_password:
            mongo_uri = mongo_uri.replace("mongodb://", f"mongodb://{mongo_user}:{mongo_password}@")
        try:
            self.mongo_client: MongoClient = MongoClient(mongo_uri)
            self.mongo_client["admin"].command("ping")
            self.mongo_db = self.mongo_client[mongo_db_name]
            self.mongo_db.command("ping")
            logger.info("Connected to MongoDB database")
            self.reaction_collection = self.mongo_db.get_collection(self.REACTION_COLLECTION_NAME)
            self.substance_collection = self.mongo_db.get_collection(self.SUBSTANCE_COLLECTION_NAME)
        except Exception as e:
            logger.error("Failed to connect to MongoDB database", exc_info=True)
            raise e

        if self.graph_backend == "memgraph":  # Memgraph
            if memgraph_user is None or memgraph_password is None:
                self._connect_memgraph(
                    uri=memgraph_uri,
                    user=memgraph_user,
                )
            else: 
                self._connect_memgraph(
                    uri=memgraph_uri,
                    user=memgraph_user,
                    password=memgraph_password,
                    encrypted=memgraph_conn_encrypted,
                )
        else:
            raise ValueError(f"Unsupported graph backend: {self.graph_backend}")

    def _connect_memgraph(
        self,
        uri: str,
        user: str | None = None,
        password: str | None = None,
        encrypted: bool = False,
    ):
        """
        Connects to a Memgraph database.
        """
        auth = None if user is None or password is None else (user, password)
        timeout = 60
        try:
            self.driver = GraphDatabase.driver(
                uri,
                auth=auth,
                max_connection_lifetime=timeout,
                encrypted=encrypted,
                trust="TRUST_ALL_CERTIFICATES",
            )
            logger.info(f"Connected to Memgraph database: {uri}")
        except Exception as e:
            logger.error(f"Failed to connect to Memgraph database: {uri}", exc_info=True)
            raise e


    def _connect_neo4j(
        self,
        uri: str,
        user: str,
        password: str,
        encrypted: bool = False,
    ):
        """
        Connects to a Neo4j database.
        """
        try:
            self.driver = GraphDatabase.driver(
                uri,
                auth=(user, password),
                encrypted=encrypted,
            )
            logger.info(f"Connected to Neo4j database: {uri}")
        except Exception as e:
            logger.error(f"Failed to connect to Neo4j database: {uri}", exc_info=True)
            raise e

    def verify_connection(self) -> bool:
        """
        Verifies the connection to the Memgraph and MongoDB databases by running simple queries.

        Returns:
            bool: True if the connection is successful, False otherwise.
        """
        return self.verify_mongo_connection() and self.verify_graphdb_connection()

    def verify_mongo_connection(self) -> bool:
        """
        Verifies the connection to the MongoDB databas by running simple queries.

        Returns:
            bool: True if the connection is successful, False otherwise.
        """
        # Verify MongoDB connection
        try:
            self.mongo_db.command("ping")
            logger.debug("MongoDB connection verified.")
            mongo_connected = True
        except Exception:
            logger.error("Failed to verify MongoDB connection.", exc_info=True)
            mongo_connected = False
        return mongo_connected

    def _verify_connection(self, session, query: str, expected_key: str, query_timeout: int) -> bool:
        """
        Helper method to execute the query and check if the result is valid.

        Args:
            session: The session object (Memgraph or Neo4j).
            query (str): The query to be executed.
            expected_key (str): The key in the result to validate.
            query_timeout (int): Timeout for the query.

        Returns:
            bool: True if the query returns the expected result, False otherwise.
        """
        try:
            logger.debug(f"Executing query: {query}")
            result = session.run(query, query_timeout=query_timeout)
            record = result.single()

            if record and expected_key in record.keys():
                return True
            else:
                logger.warning(f"No result returned or key '{expected_key}' missing.")
                return False
        except Exception as e:
            logger.error(f"Error executing query: {e}")
            return False

    def verify_graphdb_connection(self, quick_query: bool = False) -> bool:
        """
        Verifies the connection to the selected graph database (either Memgraph or Neo4j)
        based on the value of the 'graph_backend' property.

        Args:
            quick_query (bool): If True, run a lightweight query for quick validation.

        Returns:
            bool: True if the connection is successful, False otherwise.
        """
        try:
            if self.graph_backend == "memgraph":
                # Memgraph-specific connection logic
                with self.driver.session() as session:
                    query = "RETURN 1 AS result;" if quick_query else "MATCH (n) RETURN count(n) AS count"
                    expected_key = "result" if quick_query else "count"
                    query_timeout = 3 if quick_query else 30
                    return self._verify_connection(session, query, expected_key, query_timeout)

            elif self.graph_backend == "neo4j":
                # Neo4j-specific connection logic

                auth = (self.neo4j_user, self.neo4j_password) if self.neo4j_user and self.neo4j_password else None
                driver = GraphDatabase.driver(self.neo4j_uri, auth=auth, encrypted=self.neo4j_conn_encrypted)

                with driver.session() as session:
                    query = "RETURN 1 AS result;" if quick_query else "MATCH (n) RETURN count(n) AS count"
                    expected_key = "result" if quick_query else "count"
                    query_timeout = 3 if quick_query else 30
                    return self._verify_connection(session, query, expected_key, query_timeout)

            else:
                logger.warning(f"Unsupported graph backend: {self.graph_backend}")
                return False
        except Exception as e:
            logger.error(f"Error verifying graph database connection: {e}")
            return False

    def close(self):
        """
        Closes the connection to the Memgraph and MongoDB databases.
        """
        logger.debug("Closing Memgraph connection...")
        self.driver.close()
        logger.debug("Memgraph connection closed.")
        logger.debug("Closing MongoDB connection...")
        self.mongo_client.close()
        logger.debug("MongoDB connection closed.")

    @lru_cache(maxsize=1)
    def get_data_counts(self) -> Dict[str, int]:
        """
        Retrieves the total number of nodes, edges, reactions, and substances in the selected graph database (Memgraph or Neo4j) and MongoDB.

        Returns:
            Dict[str, int]: A dictionary containing the total number of nodes, edges, reactions, and substances in the databases.
        """
        return {
            "total_nodes": self._count_node_types(),
            "total_edges": self._count_edge_types(),
            "total_reactions": self._count_node_types("Reaction"),
            "total_substances": self._count_node_types("Substance"),
            "total_product_of": self._count_edge_types("PRODUCT_OF"),
            "total_reactant_of": self._count_edge_types("REACTANT_OF"),
            "total_reagent_of": self._count_edge_types("REAGENT_OF"),
            "mongo_reactions": self.reaction_collection.count_documents({}),
            "mongo_substances": self.substance_collection.count_documents({})
        }

    def _count_node_types(self, node_type: str = "") -> int:
        """
        Count the nodes in the selected graph database based on the node type.

        Args:
            node_type (str): The type of node to count.

        Returns:
            int: The count of nodes in the selected graph database.
        """
        node_type_clean = ":" + node_type if node_type else ""

        query = f"MATCH (n{node_type_clean}) RETURN count(n) AS count"
        logger.debug(f"Getting total number of {node_type} nodes in database with query '{query}'")
        with self.driver.session() as session:
            result = session.run(query)
            record = result.single()
            total_nodes = record["count"] if record else 0
            logger.debug(f"Total number of {node_type} nodes in database: {total_nodes}")

        return total_nodes

    def _count_edge_types(self, edge_type: str = "") -> int:
        """
        Count the edges in the selected graph database based on the edge type.

        Args:
            edge_type (str): The type of edge to count.

        Returns:
            int: The count of edges in the selected graph database.
        """
        edge_type_clean = ":" + edge_type if edge_type else ""

        query = f"MATCH ()-[r{edge_type_clean}]->() RETURN count(r) AS count"
        logger.debug(f"Getting total number of {edge_type} edges in database with query '{query}'")
        with self.driver.session() as session:
            result = session.run(query)
            record = result.single()
            total_edges = record["count"] if record else 0
            logger.debug(f"Total number of {edge_type} edges in database: {total_edges}")

        return total_edges

    ####################
    # Reaction methods
    ####################

    def get_reaction_by_rxid(self, rxid: str) -> Dict[str, Any]:
        """
        Retrieves a reaction node from the Mongo database by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve.

        Returns:
            Optional[Dict[str, Any]]: A dictionary containing the properties of the reaction node if found, None otherwise.

        Raises:
            ReactionNotFoundError: If the reaction with the given rxid is not found in the database.
        """
        query = {"rxid": rxid}
        reaction = self.mongo_db.get_collection(self.REACTION_COLLECTION_NAME).find_one(query, {"reagents": 0})
        if reaction:
            if reaction["smiles_is_encrypted"] in [True, "True", "true", 1]:
                if not self.SMILES_ENCRYPTION_KEY:
                    logger.error("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
                    raise ValueError("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
                try:
                    reaction["rxsmiles"] = decrypt_value_with_prepended_iv_aes_gcm(reaction["rxsmiles"], self.SMILES_ENCRYPTION_KEY)
                    # TODO : Parse extension from smiles
                    if "extension" in reaction and len(reaction["extension"]) > 0:
                        reaction["extension"] = decrypt_value_with_prepended_iv_aes_gcm(reaction["extension"], self.SMILES_ENCRYPTION_KEY)
                except Exception as e:
                    logger.error(f"Error decrypting SMILES for reaction {rxid}: {e}")
                    raise e

            return reaction
        else:
            raise ReactionNotFoundError(rxid)

    def get_reactions_paginated(self, skip: int, limit: int) -> List[Dict[str, Any]]:
        """
        Retrieves a paginated list of reaction nodes from the Mongo database.

        Args:
            skip (int): The number of reactions to skip.
            limit (int): The maximum number of reactions to retrieve.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the reaction nodes.
        """
        reactions = list(self.reaction_collection.find({}, {"reagents": 0}).skip(skip).limit(limit))
        for reaction in reactions:
            if reaction["smiles_is_encrypted"] in [True, "True", "true", 1]:
                if not self.SMILES_ENCRYPTION_KEY:
                    logger.error("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
                    raise ValueError("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
                reaction["rxsmiles"] = decrypt_value_with_prepended_iv_aes_gcm(reaction["rxsmiles"], self.SMILES_ENCRYPTION_KEY)
                # TODO : Parse extension from smiles
                # reaction["extension"] = decrypt_value_with_prepended_iv_aes_gcm(reaction["extension"], self.SMILES_ENCRYPTION_KEY)
        return reactions

    def get_reactions_count(self) -> int:
        """
        Retrieves the total number of reaction nodes in the Mongo database.

        Returns:
            int: The total number of reaction nodes.
        """
        return self.reaction_collection.count_documents({})

    def get_reactions_by_rxid(self, rxids: List[str], error_on_not_found: bool) -> Dict[str, Optional[Dict[str, Any]]]:
        """
        Retrieves reaction nodes from the Mongo database by their rxids.

        Args:
            rxids (List[str]): A list of rxids for the reactions to retrieve.

        Returns:
            Dict[str, Optional[Dict[str, Any]]]: A dictionary containing the properties of the reaction nodes for each rxid in the list,
                with the rxid as the key and the reaction properties as the value. If a reaction is not found, the value will be None.
        """
        reactions_dict: Dict[str, Optional[Dict[str, Any]]] = {}
        for rxid in rxids:
            try:
                reactions_dict[rxid] = self.get_reaction_by_rxid(rxid)
            except Exception:
                if error_on_not_found:
                    raise
                reactions_dict[rxid] = None
        return reactions_dict

    def get_rxsmiles_extended(self, rxid: str) -> str:
        """
        Retrieves the extended reaction SMILES (rxsmiles) from the Memgraph database by the rxid of the reaction.

        Args:
            rxid (str): The rxid of the reaction to retrieve the rxsmiles for.

        Returns:
            str: The extended reaction SMILES (rxsmiles) for the reaction with the given rxid.

        Raises:
            ReactionNotFoundError: If the reaction with the given rxid is not found in the database.
        """
        reaction = self.get_reaction_by_rxid(rxid)
        if "rxsmiles" not in reaction:
            raise ValueError(f"Reaction {rxid} does not have 'rxsmiles' property.")
        return reaction["rxsmiles"]

    def get_multiple_rxsmiles_extended(self, rxids: List[str]) -> Dict[str, Optional[str]]:
        """
        Retrieves the extended reaction SMILES (rxsmiles) from the Mongo database for a list of rxids.

        Args:
            rxids (List[str]): A list of rxids for the reactions to retrieve the rxsmiles for.

        Returns:
            Dict[str, Optional[str]]: A dictionary containing the rxsmiles for each reaction in the list of rxids,
                with the rxid as the key and the rxsmiles as the value. If a reaction is not found, the value will be None.
        """
        rxsmiles_dict: Dict[str, Optional[str]] = {}
        for rxid in rxids:
            try:
                rxsmiles_dict[rxid] = self.get_rxsmiles_extended(rxid)
            except Exception:
                rxsmiles_dict[rxid] = None
        return rxsmiles_dict

    def find_target_molecule_node(self, target_molecule: str) -> str:
        """
        Finds the node ID of the target molecule in the graph database by its InChIKey.

        Args:
            target_molecule (str): The InChIKey of the target molecule to find the node ID for.

        Returns:
            str: The node ID of the target molecule in the graph database.

        Raises:
            InchikeyNotFoundError: If the target molecule with the given InChIKey is not found in the database.
            MultipleInchikeyFoundError: If multiple nodes with the same InChIKey are found in the database.
        """
        query = """
                MATCH (n:Substance {inchikey: $inchikey})
                RETURN ID(n) as id
                """
        with self.driver.session() as session:
            logger.debug(f"Finding node ID for target molecule with InChIKey: {target_molecule}, query: {query}")
            result = session.run(query, inchikey=target_molecule)
            target_nodes = [record["id"] for record in result]
        if len(target_nodes) == 0:
            logger.error(f"Target molecule not found with InChIKey: {target_molecule}")
            raise InchikeyNotFoundError(target_molecule)
        elif len(target_nodes) > 1:
            logger.error(f"Multiple nodes found with InChIKey: {target_molecule}")
            raise MultipleInchikeyFoundError(target_molecule)
        else:
            logger.debug(f"Found node ID for target molecule with InChIKey: {target_molecule}")
            return str(target_nodes[0])

    def query_rxid_rxname_rxclass(self, rxid: str) -> Dict[str, str]:
        """
        Queries the knowledge base for a reaction node and its associated RXName and RXClass nodes.
        Args:
            rxid (str): The reaction ID to query.
        Returns:
            Dict[str, str]: A dictionary containing rxclass and rxname information.
        """
        query = f"""
        MATCH (c:RXClass)-[:RXCL2RXNM]->(n:RXName)-[:RXNM2R]->(r:Reaction {{rxid: '{rxid}'}})
        RETURN c, n, r
        """
        with self.driver.session() as session:
            result = session.run(query)
            record = result.single()
            if record is None:
                return {"rxname": "", "rxclass": ""}
            rxname = record.get("n", {}).get("rxname", "")
            rxclass = record.get("c", {}).get("rxclass", "")
            return {
                "rxname": rxname,
                "rxclass": rxclass,
            }

    @lru_cache(maxsize=1)
    def fetch_ambiguous_reactions(self) -> Set[str]:
        """
        Fetches the set of ambiguous reactions from the graph database.

        Returns:
            Set[str]: A set of rxids for the ambiguous reactions.
        """

        query = "MATCH (s1:Substance)-[]-(r:Reaction)-[]-(s2:Substance) WHERE s1.inchikey=s2.inchikey RETURN DISTINCT r.rxid AS rxid"

        with self.driver.session() as session:
            logger.debug(f"Fetching ambiguous reactions with query: {query}")
            results = session.run(query)
            ambig_rxns = [record["rxid"] for record in results]
            logger.debug(f"Found {len(ambig_rxns)} ambiguous reactions")

        return set(ambig_rxns)

    @lru_cache(maxsize=1)
    def fetch_multiproduct_reactions(self) -> Set[str]:
        """
        Fetches the set of multiproduct reactions from the graph database.

        Returns:
            Set[str]: A set of rxids for the multiproduct reactions.
        """

        query = "MATCH (r:Reaction)-[rel:PRODUCT_OF]->() WITH r, COUNT(rel) AS relCount WHERE relCount > 1 RETURN r.rxid AS rxid"

        with self.driver.session() as session:
            logger.debug(f"Fetching multiproduct reactions with query: {query}")
            results = session.run(query)
            multiproduct_rxns = [record["rxid"] for record in results]
            logger.debug(f"Found {len(multiproduct_rxns)} multiproduct reactions")

        return set(multiproduct_rxns)

    def get_reaction_products(self, rxid: str) -> List[Dict[str, Any]]:
        """
        Retrieves the products of a reaction node from the graph database by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve the products for.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the product nodes for the reaction.
        """
        # Get product edges from reaction
        query = """
                MATCH (r:Reaction {rxid: $rxid})-[:PRODUCT_OF]->(s:Substance)
                RETURN ID(s) AS node_id, s.inchikey AS inchikey
                """
        with self.driver.session() as session:
            logger.debug(f"Fetching products of reaction with rxid: {rxid}, query: {query}")
            result = session.run(query, rxid=rxid)
            products_nodes = [record for record in result]
            logger.debug(f"Found {len(products_nodes)} products for reaction with rxid: {rxid}")

        # Search for each substance by inchikey in mongo
        products = []
        for node in products_nodes:
            inchikey = node["inchikey"]
            substance = self.get_substance_by_inchikey(inchikey)
            # Manually add node_id to substance
            # TODO - Ideally this mapping could already exist in MongoDB somehow?
            substance["graph_node_id"] = node["node_id"]
            products.append(substance)

        return products

    def get_reaction_reactants(self, rxid: str) -> List[Dict[str, Any]]:
        """
        Retrieves the reactants of a reaction node from the graph database by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve the reactants for.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the reactant nodes for the reaction.
        """
        # Get reactant edges from reaction
        query = """
                MATCH (r:Reaction {rxid: $rxid})<-[:REACTANT_OF]-(s:Substance)
                RETURN ID(s) AS node_id, s.inchikey AS inchikey
                """
        with self.driver.session() as session:
            logger.debug(f"Fetching reactants of reaction with rxid: {rxid}, query: {query}")
            result = session.run(query, rxid=rxid)
            reactants_nodes = [record for record in result]
            logger.debug(f"Found {len(reactants_nodes)} reactants for reaction with rxid: {rxid}")

        # Search for each substance by inchikey in mongo
        reactants = []
        for node in reactants_nodes:
            inchikey = node["inchikey"]
            substance = self.get_substance_by_inchikey(inchikey)
            # Manually add node_id to substance
            # TODO - Ideally this mapping could already exist in MongoDB somehow?
            substance["graph_node_id"] = node["node_id"]
            reactants.append(substance)

        return reactants

    def get_reaction_reagents(self, rxid: str) -> List[Dict[str, Any]]:
        """
        Retrieves the reagents of a reaction node from the graph database by its rxid.

        Args:
            rxid (str): The rxid of the reaction to retrieve the reagents for.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the reagent nodes for the reaction.
        """
        # Get reagent edges from reaction
        query = """
                MATCH (r:Reaction {rxid: $rxid})<-[:REAGENT_OF]-(s:Substance)
                RETURN ID(s) AS node_id, s.inchikey AS inchikey
                """
        with self.driver.session() as session:
            logger.debug(f"Fetching reagents of reaction with rxid: {rxid}, query: {query}")
            result = session.run(query, rxid=rxid)
            reagents_nodes = [record for record in result]
            logger.debug(f"Found {len(reagents_nodes)} reagents for reaction with rxid: {rxid}")

        # Search for each substance by inchikey in mongo
        reagents = []
        for node in reagents_nodes:
            inchikey = node["inchikey"]
            substance = self.get_substance_by_inchikey(inchikey)
            # Manually add node_id to substance
            # TODO - Ideally this mapping could already exist in MongoDB somehow?
            substance["graph_node_id"] = node["node_id"]
            reagents.append(substance)

        return reagents

    #####################
    # Substance methods
    #####################

    def get_substances_paginated(self, skip: int, limit: int) -> List[Dict[str, Any]]:
        """
        Retrieves a paginated list of substance nodes from the Mongo database.

        Args:
            skip (int): The number of substances to skip.
            limit (int): The maximum number of substances to retrieve.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the properties of the substance nodes.
        """
        substances = list(self.substance_collection.find({}, {"reagents": 0}).skip(skip).limit(limit))
        for substance in substances:
            # TODO : Add smiles_is_encrypted to substance
            # if substance["smiles_is_encrypted"] in [True, "True", "true", 1]:
            if not self.SMILES_ENCRYPTION_KEY:
                logger.error("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
                raise ValueError("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
            try:
                substance["canonical_smiles"] = decrypt_value_with_prepended_iv_aes_gcm(substance["canonical_smiles"], self.SMILES_ENCRYPTION_KEY)
            except Exception as e:
                logger.error(f"Error decrypting SMILES for substance {substance['inchikey']}: {e}")
                raise e
        return substances

    def get_substances_count(self) -> int:
        """
        Retrieves the total number of substance nodes in the Mongo database.

        Returns:
            int: The total number of substance nodes.
        """
        return self.substance_collection.count_documents({})

    def get_substance_by_inchikey(self, inchikey: str) -> Dict[str, Any]:
        """
        Retrieves a substance node from the Mongo database by its InChIKey.

        Args:
            inchikey (str): The InChIKey of the substance to retrieve.

        Returns:
            Optional[Dict[str, Any]]: A dictionary containing the properties of the substance node if found, None otherwise.
        """
        query = {"inchikey": inchikey}
        substance = self.substance_collection.find_one(query, {"reagents": 0})
        if substance is not None:
            # TODO : Add smiles_is_encrypted to substance
            # if substance["smiles_is_encrypted"] in [True, "True", "true", 1]:
            if not self.SMILES_ENCRYPTION_KEY:
                logger.error("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
                raise ValueError("SMILES_ENCRYPTION_KEY is not set in the environment. Cannot decrypt SMILES.")
            try:
                substance["canonical_smiles"] = decrypt_value_with_prepended_iv_aes_gcm(substance["canonical_smiles"], self.SMILES_ENCRYPTION_KEY)
            except Exception as e:
                logger.error(f"Error decrypting SMILES for substance {substance['inchikey']}: {e}")
                raise e
            return substance
        raise InchikeyNotFoundError(inchikey)

    def get_substance_product_in(self, inchikey: str) -> List[str]:
        """
        Retrieves a list of reaction IDs where the substance with the given InChIKey is a product.

        Args:
            inchikey (str): The InChIKey of the substance to retrieve the reactions for.

        Returns:
            List[str]: A list of reaction IDs where the substance with the given InChIKey is a product.
        """
        query = """
                MATCH (s:Substance {inchikey: $inchikey})<-[:PRODUCT_OF]-(r:Reaction)
                RETURN r.rxid AS rxid
                """

        with self.driver.session() as session:
            result = session.run(query, inchikey=inchikey)
            reaction_nodes = [record for record in result]
        reactions = []
        for node in reaction_nodes:
            rxid = node["rxid"]
            reactions.append(rxid)
        return reactions

    def get_substance_reactant_in(self, inchikey: str) -> List[str]:
        """
        Retrieves a list of reaction IDs where the substance with the given InChIKey is a reactant.

        Args:
            inchikey (str): The InChIKey of the substance to retrieve the reactions for.

        Returns:
            List[str]: A list of reaction IDs where the substance with the given InChIKey is a reactant.
        """
        query = """
                MATCH (s:Substance {inchikey: $inchikey})-[:REACTANT_OF]->(r:Reaction)
                RETURN r.rxid AS rxid
                """

        with self.driver.session() as session:
            result = session.run(query, inchikey=inchikey)
            reaction_nodes = [record for record in result]
        reactions = []
        for node in reaction_nodes:
            rxid = node["rxid"]
            reactions.append(rxid)
        return reactions

    def get_substance_reagent_in(self, inchikey: str) -> List[str]:
        """
        Retrieves a list of reaction IDs where the substance with the given InChIKey is a reagent.

        Args:
            inchikey (str): The InChIKey of the substance to retrieve the reactions for.

        Returns:
            List[str]: A list of reaction IDs where the substance with the given InChIKey is a reagent.
        """
        query = """
                MATCH (s:Substance {inchikey: $inchikey})-[:REAGENT_OF]->(r:Reaction)
                RETURN r.rxid AS rxid
                """

        with self.driver.session() as session:
            result = session.run(query, inchikey=inchikey)
            reaction_nodes = [record for record in result]
        reactions = []
        for node in reaction_nodes:
            rxid = node["rxid"]
            reactions.append(rxid)
        return reactions

    #################
    # Graph methods
    #################
    def fetch_synthesis_graph(self, search_params: SynthGraphSearch) -> SynthGraph:
        """
        Fetches a synthesis graph from the graph database based on the given search parameters.

        Args:
            search_params (SynthGraphSearch): The search parameters to use for fetching the synthesis graph.

        Returns:
            SynthGraph: A SynthGraph object containing the synthesis graph and other relevant information.

        Raises:
            ValueError: If the search depth is not an integer.
            InchikeyNotFoundError: If the target molecule with the given InChIKey is not found in the database.
        """
        query_type = search_params.query_type
        # Ensure 'depth' is an integer to prevent injection
        if not isinstance(search_params.search_depth, int):
            raise ValueError("Depth must be an integer")

        if query_type == "shortest_path":
            if search_params.graph_backend == "neo4j":
                query = f"""
                    MATCH (ra:Substance)-[:PRODUCT_OF|REAGENT_OF|REACTANT_OF*..{search_params.search_depth}]->(pr:Substance {{inchikey: $target_molecule}})
                    WHERE ra.inchikey <> pr.inchikey
                    WITH COLLECT(DISTINCT ra) AS sms, pr
                    UNWIND sms AS sm
                    MATCH p=shortestPath((sm)-[:PRODUCT_OF|REAGENT_OF|REACTANT_OF*..{search_params.search_depth}]->(pr))
                    WITH [n IN nodes(p) WHERE n:Reaction] AS rxns
                    UNWIND rxns AS rxn
                    MATCH (rxn)-[rel:REACTANT_OF|REAGENT_OF|PRODUCT_OF*1..1]-(s_final:Substance)
                    RETURN rxn, s_final, rel
                """
            elif search_params.graph_backend == "memgraph":
                query = (
                    f"MATCH (ra:Substance)-[frels:PRODUCT_OF|REACTANT_OF*..{search_params.search_depth}]->(pr:Substance {{inchikey: $target_molecule}}) "
                    "WHERE ra.inchikey <> pr.inchikey "
                    "WITH COLLECT(DISTINCT ra) AS sms, pr "
                    "UNWIND sms AS sm "
                    f"MATCH p = (sm)-[relationships:PRODUCT_OF|REACTANT_OF *BFS ..{search_params.search_depth}]->(pr) "
                    "UNWIND nodes(p) AS rxn "
                    "MATCH (rxn:Reaction)-[rel:REACTANT_OF|REAGENT_OF|PRODUCT_OF]-(s_final:Substance) "
                    "RETURN rxn, s_final, rel;"
                )
        elif query_type == "full_graph":
            query = f"""
                MATCH p = (ra:Substance)-[:PRODUCT_OF|REAGENT_OF|REACTANT_OF*..{search_params.search_depth}]->(pr:Substance {{inchikey: $target_molecule}})
                WHERE pr.inchikey <> ra.inchikey RETURN p
            """
        else:
            # TODO : Create more specific error
            raise Exception("Unsupported query type")

        G = DiGraph()

        with self.driver.session() as session:
            try:
                logger.debug(f"Fetching synthesis graph with query: {query}")
                resulting_graph = session.run(query, target_molecule=search_params.target_molecule_inchikey).graph()
            except Exception:
                logger.error("Failed to fetch synthesis graph", exc_info=True)
                raise Exception("Failed to fetch synthesis graph")

            logger.debug(
                f"Synthesis graph fetched with {len(resulting_graph.nodes)} nodes and {len(resulting_graph.relationships)} edges for target molecule: {search_params.target_molecule_inchikey}"
            )
            for node in resulting_graph.nodes:
                node_properties = self._extract_node_properties(node)
                G.add_node(node.id, **node_properties)

            for relationship in resulting_graph.relationships:
                edge_properties = self._extract_edge_properties(relationship)
                G.add_edge(
                    relationship.start_node.id,
                    relationship.end_node.id,
                    **edge_properties,
                )

        # Find target node in networkx graph
        target_node = None
        for node_id in G.nodes:
            # TODO : Check if this is needed for networkx graphs
            # Set node_id for each node in the graph (if necessary)
            G.nodes[node_id]["node_id"] = node_id

            # Process substances
            if G.nodes[node_id].get("node_type") == "substance":
                inchikey = G.nodes[node_id]["inchikey"]
                if inchikey == search_params.target_molecule_inchikey:
                    target_node = node_id

        if target_node is None:
            logger.error(f"Target molecule not found in synthesis graph: {search_params.target_molecule_inchikey}")
            raise SubstanceNotFoundInSynthGraphError(search_params.target_molecule_inchikey)

        synth_graph = SynthGraph(
            target_molecule_node_id=target_node,
            synthesis_graph=G.copy(),
            search_params=search_params,
        )

        return synth_graph

    def _extract_node_properties(self, node) -> Dict[str, Any]:
        """
        Extracts the properties of a node from the graph database.

        Args:
            node: The node to extract the properties from.

        Returns:
            Dict[str, Any]: A dictionary containing the properties of the node.

        Raises:
            Exception: If the node type is not recognized.
        """
        node_properties = {}
        for k, v in node.items():
            node_properties[k] = v

        if node_properties.get("node_type", None) is None:
            label = list(node.labels)[0] or "unknown"
            node_properties["node_type"] = label.lower()

        node_type = node_properties.get("node_type", None).lower()
        if node_type == "substance":
            node_properties["node_label"] = node_properties["inchikey"]
            if "canonical_smiles" in node_properties:
                node_properties["canonical_smiles"] = decrypt_value_with_prepended_iv_aes_gcm(
                    node_properties.get("canonical_smiles", ""), self.SMILES_ENCRYPTION_KEY
                )
        elif node_type == "reaction":
            node_properties["node_label"] = node_properties["rxid"]
            if "rxsmiles" in node_properties:
                node_properties["rxsmiles"] = decrypt_value_with_prepended_iv_aes_gcm(node_properties.get("rxsmiles", ""), self.SMILES_ENCRYPTION_KEY)
            node_properties["yield_info"] = {
                "yield_predicted": node_properties.get("yield_predicted", 0.0),
                "yield_score": node_properties.get("yield_score", 0.0),
            }
        else:
            logger.debug(f"Node type not recognized: {node_type}")
            raise Exception("[ERROR] Invalid node type.")
        return node_properties

    def _extract_edge_properties(self, relationship) -> Dict[str, Any]:
        """
        Extracts the properties of an edge from the graph database.

        Args:
            relationship: The relationship to extract the properties from.

        Returns:
            Dict[str, Any]: A dictionary containing the properties of the edge.

        Raises:
            Exception: If the edge type is not recognized.
        """
        edge_properties = {}
        for k, v in relationship.items():
            edge_properties[k] = v

        if edge_properties.get("edge_type", None) is None:
            label = relationship.type or "unknown"
            edge_properties["edge_type"] = label.lower()

        forward_s2r_edge_types = ["reactant_of", "reagent_of"]
        reverse_s2r_edge_types = ["product_of"]

        start_node = relationship.start_node
        end_node = relationship.end_node
        edge_properties["start_node_id"] = start_node.id
        edge_properties["end_node_id"] = end_node.id
        if edge_properties.get("edge_type") in forward_s2r_edge_types:
            edge_properties["inchikey"] = start_node["inchikey"]
            edge_properties["rxid"] = end_node["rxid"]
            edge_properties["start_node"] = start_node["inchikey"]
            edge_properties["end_node"] = end_node["rxid"]
        elif edge_properties.get("edge_type") in reverse_s2r_edge_types:
            edge_properties["inchikey"] = end_node["inchikey"]
            edge_properties["rxid"] = start_node["rxid"]
            edge_properties["start_node"] = start_node["rxid"]
            edge_properties["end_node"] = end_node["inchikey"]
        else:
            # TODO : Create more specific error
            raise Exception("[ERROR] Invalid edge type.")

        return edge_properties
