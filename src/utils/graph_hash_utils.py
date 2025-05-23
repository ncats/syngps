import hashlib
import json
from networkx import DiGraph
import logging

logger = logging.getLogger(__name__)

def hash_graph(graph: DiGraph) -> str:
    """
    Create a hash for a graph, ensuring all node and edge IDs are treated as strings.

    Args:
        graph (DiGraph): The graph to hash.

    Returns:
        str: The hash of the graph.
    """
    # Extract node and edge IDs, convert them to strings, and sort
    node_ids = sorted(str(node) for node in graph.nodes())
    edge_ids = sorted((str(u), str(v)) for u, v in graph.edges())

    # Serialize the structure with nodes and edges as strings
    graph_structure = json.dumps({"nodes": node_ids, "edges": edge_ids}, sort_keys=True)

    logger.debug(f"Hashing graph structure: {graph_structure}")
    hash = hashlib.sha256(graph_structure.encode()).hexdigest()
    logger.debug(f"Graph hash: {hash}")

    # Create a hash of the serialized structure
    return hash