# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: NCATS/NIH


# References



import requests
import json
import pandas as pd
import networkx as nx


# Definitions


# Cytoscape API URL
CYTOSCAPE_API_URL = "http://localhost:1234/v1"

# new_style_name = "default"





LAYOUT_TYPE = 'hierarchical'

NEW_STYLE_NAME = 'BOYA'

NEW_STYLE_JSON = {
    "title": NEW_STYLE_NAME,
    "defaults": [
        {"visualProperty": "NODE_SIZE", "value": 40},
        {"visualProperty": "EDGE_LINE_TYPE", "value": "SOLID"},
        {"visualProperty": "EDGE_WIDTH", "value": 2},
        {"visualProperty": "EDGE_CURVED", "value": False},
        {"visualProperty": "EDGE_TARGET_ARROW_SHAPE", "value": "DELTA"}
        
    ],
    "mappings": [
        {
            "mappingType": "discrete",
            "mappingColumn": "node_type",
            "mappingColumnType": "String",
            "visualProperty": "NODE_FILL_COLOR",
            "map": [
                {"key": "reaction", "value": "#F19983"}
                
            ]
            
        },
        {
            "mappingType": "discrete",
            "mappingColumn": "node_type",
            "mappingColumnType": "String",
            "visualProperty": "NODE_SHAPE",
            "map": [
  
                {"key": "reaction", "value": "ELLIPSE"}

                
            ]
            
        },   
        {
            "mappingType": "discrete",
            "mappingColumn": "edge_type",
            "mappingColumnType": "String",
            "visualProperty": "EDGE_STROKE_UNSELECTED_PAINT",
            "map": [
                {"key": "reverse_yield", "value": "#000000"}
                
            ]
            
        },


        {
            "mappingType": "passthrough",
            "mappingColumn": "node_id",
            "mappingColumnType": "String",
            "visualProperty": "NODE_LABEL"
            
        }

    ]
}






def route2cyjs (G):

    cy_json_data = nx.cytoscape_data(G, name='uuid', ident='node_label')

    return (cy_json_data)


def show_in_cytotscape(cy_json):
    
        
    # Send the network to Cytoscape
    response = requests.post(f"{CYTOSCAPE_API_URL}/networks?format=cyjs", json = cy_json)
    SUID = None
    
    
    if response.ok:
        network_suid = response.json()['networkSUID']
        print(f"Network created with SUID: {network_suid}")
        
        # Create a view for the network
        view_response = requests.get(f"{CYTOSCAPE_API_URL}/networks/{network_suid}/views/first")
        if view_response.ok:
            print("Network view created.")
            
            SUID = int(view_response.json()['data']['SUID'])

            return (SUID)
        else:
            print("Failed to create network view.")

        return (None)
    else:
        print("Failed to create network.")

    return (None)



def create_style ():

    
    # Create style if does not exist
    
    

    

    
    # Check if the style already exists
    existing_styles_response = requests.get(f"{CYTOSCAPE_API_URL}/styles")
    if existing_styles_response.ok:
        existing_styles = existing_styles_response.json()
        print("Existing styles:", existing_styles)  # Debug print
    
        style_names = [style['title'] for style in existing_styles if isinstance(style, dict)]
        
        if NEW_STYLE_NAME in style_names:
            print(f"Style '{NEW_STYLE_NAME}' already exists. Applying existing style.")
        else:
            print(f"Creating new style '{NEW_STYLE_NAME}'.")
    
    
    
            # Create the new style in Cytoscape
            create_style_response = requests.post(f"{CYTOSCAPE_API_URL}/styles", json = NEW_STYLE_JSON)
    
            if create_style_response.ok:
                print(f"New style '{NEW_STYLE_NAME}' created.")

                return (True)
            
            else:
                print("Failed to create new style.")
    else:
        print("Failed to retrieve existing styles.")

    return (False)


def apply_style (network_suid):

    create_style ()

    
    
    apply_style_response = requests.get(f"{CYTOSCAPE_API_URL}/apply/styles/{NEW_STYLE_NAME}/{network_suid}")
    #print (apply_style_response)
    
    if apply_style_response.ok:
        print(f"Style '{NEW_STYLE_NAME}' applied to the network.")
        
        return (True)
    else:
        print(f"Failed to apply style '{NEW_STYLE_NAME}'.")

    return (False)






def apply_layout (network_suid):

    apply_style_response = requests.get(f"{CYTOSCAPE_API_URL}/apply/layouts/{LAYOUT_TYPE}/{network_suid}")
    
    if apply_style_response.ok:
        print(f"Layout '{LAYOUT_TYPE}' applied to the network.")
        
        return (True)
    else:
        print(f"Failed to apply layout '{LAYOUT_TYPE}'.")

    return (False)    


def rename_network (nw_name, suuid):



    
    
    #rename_response = requests.post(f"{CYTOSCAPE_API_URL}/commands/network/rename", json = {"name": nw_name, "sourceNetwork": "current"})
    rename_response = requests.post(f"{CYTOSCAPE_API_URL}/commands/network/rename", json = {"name": nw_name, "sourceNetwork": f"SUID:{suuid}"})
    #print (apply_style_response)
    
    if rename_response.ok:
        print(f"Style '{nw_name}' applied to the network (SUID: {suuid}).")
        
        return (True)
    else:
        print(f"Failed to rename network (SUID: {suuid}) to name: '{nw_name}'.")

    return (False)