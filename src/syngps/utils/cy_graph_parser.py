# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: NCATS/NIH

# References
#
# Ref: https://networkx.org/documentation/stable/reference/generated/networkx.classes.function.induced_subgraph.html
#
import networkx as nx
import json
import logging
logger = logging.getLogger(__name__)

def parse_input(fname):
    """
    Just parses in the input .cyson coming Cytoscape.
    """

    with open(fname, 'r') as file:
        data = json.load(file)

    return (data)


def add_node(G, node_id, node_metadata):
    """
    Helper function to add node with its metadata to the directed NetworkX object (DiGraph).
    """
    if node_id not in G:
        G.add_node(node_id)
    else:
        logger.error(
            f'[ERROR] Node of duplicate ID: {node_id} present in network.')
        raise ValueError(f'[ERROR] Node of duplicate ID: {node_id} present in network.')

    for key in node_metadata.keys():
        G.nodes[node_id][key] = node_metadata[key]

    G.nodes[node_id]['node_id'] = node_id

    return (G)


def process_nodes(G, nodes):
    """
        Helper function to create the nodes and add metadata for the NetworkX object.
        Input:  - an (empty) NetworkX directed graph object (DiGraph),
                - nodes section of the ASKCOS JSON file,

        Output: Directed NetworkX DiGraph object.
    """
    valid_node_types = ['reaction', 'substance']

    for n in nodes:

        node_metadata = {}

        node_type = n['data']['node_type']

        if node_type not in valid_node_types:
            logger.error(
                f'[ERROR] Invalid node type found: {node_type}. Terminating ...'
            )
            raise ValueError(
                f'[ERROR] Invalid node type found: {node_type}. Terminating ...'
            )
        node_uuid = str(n['data']['SUID'])

        # Ignoring source info. Filtering should be done before this stage,
        # because at this stage we're only interested in synthesis routes.

        node_metadata['uuid'] = node_uuid

        if node_type == 'reaction':
            node_metadata['node_type'] = node_type

            rxid = n['data']['name']
            yield_score = 1.00
            yield_predicted = 1.00

            rxsmiles = 'rxsmiles_placeholder'

            node_metadata['rxid'] = rxid

            node_metadata['node_label'] = n['data']['name']
            node_metadata['name'] = n['data']['name']
            # node_metadata['uuid'] = node_uuid
            node_metadata['yield_score'] = round(yield_score, 2)
            node_metadata['yield_predicted'] = round(yield_predicted, 2)
            node_metadata['rxsmiles'] = rxsmiles
            # node_metadata['aggregated_yield'] = n['aggregated_yield']

            G = add_node(G, node_uuid, node_metadata)

        else:
            node_metadata['node_type'] = 'substance'

            smiles = 'smiles_placeholder'
            inchikey = 'inchikey_placeholder'

            node_metadata['inchikey'] = inchikey
            node_metadata['node_label'] = n['data']['name']
            node_metadata['name'] = n['data']['name']
            # node_metadata['uuid'] = n['data']['SUID']
            node_metadata['srole'] = n['data']['srole']

            # !!!for testing only!!!

            # node_metadata['srole'] = 'im'

            # placeholder for parsing synthesis planning roles of substances
            # i.e.: starting material, intermedier, target molecule

            G = add_node(G, node_uuid, node_metadata)

    return (G)


def add_edge(G, start_node, end_node, edge_metadata):
    """
        Helper function to register edges.

        Input: start node, end node, edge metadata.
        Otput: directed NetworkX object (DiGraph)
    """
    G.add_edge(start_node, end_node)

    for key in edge_metadata.keys():
        G[start_node][end_node][key] = edge_metadata[key]

    return (G)


def process_edges(G, edges):
    """
        Parsing edges.

        Input: Directed NetworkX object (DiGraph),

        Output: Synthesis graph.
    """

    for e in edges:
        edge_metadata = {}

        start_node = str(e['data']['source'])
        end_node = str(e['data']['target'])

        edge_metadata['uuid'] = str(e['data']['SUID'])
        edge_metadata['edge_type'] = e['data']['edge_type']

        G = add_edge(G, start_node, end_node, edge_metadata)

    return (G)


# By ChatGPT
# Convert the graph to JSON
def graph_to_json(graph):
    # Convert nodes
    nodes = [
        # Include metadata as key-value pairs
        {"id": node, **graph.nodes[node]}
        for node in graph.nodes
    ]
    # Convert edges
    edges = [
        {"start_node": u, "end_node": v, **
            graph.edges[u, v]}  # Include metadata
        for u, v in graph.edges
    ]
    # Combine into a single dictionary
    graph_json = {"nodes": nodes, "edges": edges}

    return (graph_json)


def add_aggregated_yield_info(G, agg_yield):
    for n in G.nodes:
        G.nodes[n]['aggregated_yield'] = agg_yield

    return (G)


def cy_json2synthgraph(data_json):

    nodes = data_json['elements']['nodes']
    edges = data_json['elements']['edges']

    G = nx.DiGraph()

    G = process_nodes(G, nodes)
    G = process_edges(G, edges)

    return (G)
