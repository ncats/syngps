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


def parse_input(fname):
    """
        Just parses in the input json coming from AICP server_url/aicp-synthplanning/api/v1/prediction/synthesis_routes EP, saved as file.
    """
    
    with open(fname, 'r') as file:
        data = json.load(file)

    return (data)





def add_node (G, node_id, node_metadata):
    """
        Helper function to add node with its metadata to the directed NetworkX object (DiGraph).
    """


    if node_id not in G:
        G.add_node(node_id)
    else:
        print (f'[ERROR] Node of duplicate ID: {node_id} present in network. Terminating...')
        sys.exit (-1)
    
    
    for key in node_metadata.keys():
        G.nodes[node_id][key] = node_metadata[key]
    
    G.nodes[node_id]['node_id'] = node_id
    
    return (G)


def process_nodes (G, nodes):
    """
        Helper function to create the nodes and add metadata for the NetworkX object.
        Input:  - an (empty) NetworkX directed graph object (DiGraph),
                - nodes section of the ASKCOS JSON file,
                
        Output: Directed NetworkX DiGraph object.


    """
    valid_node_types = ['reaction', 'substance']



    for n in nodes:

        node_metadata = {}

        node_type = n['node_type']

        if node_type not in valid_node_types:
            write_to_log ("Invalid node type found: %s. Terminating ..." % (node_type))
            sys.exit(-1)

        node_uuid = n['uuid']

        # Ignoring source info. Filtering should be done before this stage,
        # because at this stage we're only interested in synthesis routes.


        node_metadata['uuid'] = node_uuid



        if node_type == 'reaction':
            node_metadata['node_type'] = node_type

            rxid = n['rxid']
            yield_score = float(n['yield_info']['yield_score'])
            yield_predicted = float(n['yield_info']['yield_predicted'])


            rxsmiles = n['rxsmiles']


            node_metadata['rxid'] = rxid

            node_metadata['node_label'] = n['node_label']
            node_metadata['name'] = n['node_label']
            node_metadata['uuid'] = n['uuid']
            node_metadata['yield_score'] = round(yield_score, 2)
            node_metadata['yield_predicted'] = round(yield_predicted, 2)
            node_metadata['rxsmiles'] = rxsmiles
            #node_metadata['aggregated_yield'] = n['aggregated_yield']




            



            G = add_node (G, rxid, node_metadata)


            


        else:
            node_metadata['node_type'] = 'substance'


            smiles = n['canonical_smiles']
            inchikey = n['inchikey']

            node_metadata['uuid'] = n['uuid']
           


            
            node_metadata['inchikey'] = inchikey
            node_metadata['node_label'] = n['node_label']
            node_metadata['name'] = n['node_label']
            node_metadata['uuid'] = n['uuid']
            node_metadata['srole'] = n['srole']


            # !!!for testing only!!!

            #node_metadata['srole'] = 'im'    

            # placeholder for parsing synthesis planning roles of substances
            # i.e.: starting material, intermedier, target molecule

            G = add_node (G, inchikey, node_metadata)







    return (G)



def add_edge (G, start_node, end_node, edge_metadata):
    """
        Helper function to register edges.
        
        Input: start node, end node, edge metadata.
        Otput: directed NetworkX object (DiGraph)
    """
    G.add_edge (start_node, end_node)
    
    for key in edge_metadata.keys():
        G[start_node][end_node][key] = edge_metadata[key]


    return (G)




def process_edges (G, edges):
    """
        Parsing edges.

        Input: Directed NetworkX object (DiGraph),
        
        Output: Synthesis graph.


                
    """

    for e in edges:
        
        edge_metadata ={}

        #edge_uuid = e['uuid']
        #edge_type = e['edge_type']

        #if edge_type not in valid_edge_types:
        #    write_to_log ("Invalid edge type found: %s. Terminating ..." % (edge_type))
        #    sys.exit(-1)

        start_node = e['start_node']
        end_node = e['end_node']



            
        
        #edge_metadata['start_node'] = start_node
        #edge_metadata['end_node'] = end_node

        edge_metadata['node_label'] = e['edge_label']
        edge_metadata['name'] = e['edge_label']
        edge_metadata['uuid'] = e['uuid']
        edge_metadata['edge_type'] = e['edge_type']

        

        G = add_edge (G, start_node, end_node, edge_metadata)
            

    
    return (G)





# By ChatGPT
# Convert the graph to JSON
def graph_to_json(graph):
    # Convert nodes
    nodes = [
        {"id": node, **graph.nodes[node]}  # Include metadata as key-value pairs
        for node in graph.nodes
    ]
    # Convert edges
    edges = [
        {"start_node": u, "end_node": v, **graph.edges[u, v]}  # Include metadata
        for u, v in graph.edges
    ]
    # Combine into a single dictionary
    graph_json = {"nodes": nodes, "edges": edges}

    return (graph_json)


def add_aggregated_yield_info (G, agg_yield):
    for n in G.nodes:
        print (agg_yield)
        G.nodes[n]['aggregated_yield'] = agg_yield

    return (G)

def json2routes (data_json):
    # First prarse synthesis graph
    # Then create the node induced subgraphs as defined by the routes section



    
    nodes = data_json['evidence_synth_graph']['nodes']
    edges = data_json['evidence_synth_graph']['edges']
    
    synthesis_routes = []

    routes = data_json['routes']


    G = nx.DiGraph()

    
    G = process_nodes (G, nodes)
    G = process_edges (G, edges)

    #print (json.dumps(graph_to_json(G), indent = 4))

    for route in routes:
        S = None
        if route['predicted'] ==  False:
            node_ids = route['route_node_labels']
            #print (node_ids)

            H = G.copy()
            S = nx.induced_subgraph(H, node_ids)
            S = add_aggregated_yield_info (S, route['aggregate_yield'])
            
            synthesis_routes.append(S)

        #print (graph_to_json(S))

    # Adding synthesis graph as last "route"
    synthesis_routes.append(G)

    
    return (synthesis_routes)











