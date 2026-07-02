from typing import Any, Dict, Literal

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from app.adapters import AskcosApiAdapter
from syngps.models.models import SynthGraph
from app.app_config import APP_CONFIG
from app.models import (
    AtommapRequest,
    AtommapResponse,
    ModuleStatusesModel,
    TreeSearchInput,
    TreeSearchResponse,
    AskcosConditionPredictRequest,
    AskcosConditionPredictResponse,
    ComputeBalanceIndicesResponse
)
from app.utils import (
    askcos_tree2synth_graph,
    askcos_tree2synth_paths_with_graph,
    is_reaction_smiles_parseable,
    uds_tree2synth_graph,
    uds_tree2synth_paths,
    uds_tree2synth_paths_with_graph
)
from syngps import decomposition_utils
from app.logging_config import logger

askcos_adapter = AskcosApiAdapter(askcos_base_url=APP_CONFIG.askcos_base_url)

###########################
# ASKCOS Wrappers Router
###########################
askcos_wrapper_router = APIRouter(prefix="/hypothesis-engine/askcos", tags=["ASCKOS"])


@askcos_wrapper_router.get("/status", summary="Get the status of the ASKCOS service")
def get_askcos_status() -> ModuleStatusesModel:
    """
    Get the status of the ASKCOS service.

    Returns:
        dict: A dictionary with the status of the ASKCOS service
    """
    return askcos_adapter.get_askcos_status()


@askcos_wrapper_router.post(
    "/atom-map/rxnmapper", summary="Get atom mapping for a reaction using the ASKCOS RXNMapper endpoint (/api/atom-map/rxnmapper/call-sync)"
)
def get_atom_mapping_rxnmapper(request: AtommapRequest) -> AtommapResponse:
    """
    Get atom mapping for a reaction using the ASKCOS RXNMapper endpoint (/api/atom-map/rxnmapper/call-sync).

    Args:
        request (AtommapRequest): Request object containing the reaction to map

    Returns:
        AtommapResponse: Response object containing the mapped reaction
    """
    if ">" not in request.smiles:
        raise HTTPException(status_code=420, detail="Input must be reaction smiles.")

    if not is_reaction_smiles_parseable(request.smiles):
        raise HTTPException(status_code=420, detail="Reaction smiles not parseable by RDKIT.")

    return askcos_adapter.askcos_atommap_rxnmapper(request)


# @askcos_wrapper_router.post(
#     "/tree-search/raw", summary="Get raw ASKCOS tree search results, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token')"
# )
# def get_tree_search_raw(input: TreeSearchInput) -> TreeSearchResponse:
#     """
#     Get raw ASKCOS tree search results, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token').

#     Args:
#         input (dict): Input object containing the query reaction

#     Returns:
#         dict: Response object containing the raw ASKCOS response
#     """
#     return askcos_adapter.askcos_tree_search_raw(input)

@askcos_wrapper_router.post(
    "/tree-search/raw", summary="Get raw ASKCOS tree search results, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token')"
)
def get_tree_search_raw(input: TreeSearchInput) -> Any:
    """
    Get raw ASKCOS tree search results, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token').

    Args:
        input (dict): Input object containing the query reaction

    Returns:
        dict: Response object containing the raw ASKCOS response
    """
    return askcos_adapter.askcos_tree_search_raw_v2(input)


# @askcos_wrapper_router.post(
#     "/tree-search/parse-graph",
#     summary="Get ASKCOS tree search results and parses the graph in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token')",
# )
# def get_tree_search_parse_graph(input: TreeSearchInput) -> Any:
#     """
#     Get ASKCOS tree search results and parses the graph in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token').

#     Args:
#         input (dict): Input object containing the query reaction

#     Returns:
#         dict: Response object containing the parsed ASKCOS graph
#     """
#     raw_response = askcos_adapter.askcos_tree_search_raw(input)
#     response = askcos_tree2synth_graph(raw_response)
#     return SynthGraph(synthesis_graph=response, target_molecule_node_id="BLANK")

@askcos_wrapper_router.post(
    "/tree-search/parse-graph",
    summary="Get ASKCOS tree search results and parses the graph in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token')",
)
def get_tree_search_parse_graph(input: TreeSearchInput) -> Any:
    """
    Get ASKCOS tree search results and parses the graph in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token').

    Args:
        input (dict): Input object containing the query reaction

    Returns:
        dict: Response object containing the parsed ASKCOS graph
    """
    raw_response = askcos_adapter.askcos_tree_search_raw_v2(input)
    uds = raw_response["result"]["uds"]
    response = uds_tree2synth_graph(uds)
    return SynthGraph(synthesis_graph=response, target_molecule_node_id="BLANK")


# @askcos_wrapper_router.post(
#     "/tree-search/parse-paths",
#     summary="Get ASKCOS tree search results and parses the routes in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token')",
# )
# def get_tree_search_parse_paths(input: TreeSearchInput) -> Any:
#     """
#     Get ASKCOS tree search results and parses the paths in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token').

#     Args:
#         input (dict): Input object containing the query reaction

#     Returns:
#         dict: Response object containing the parsed ASKCOS paths
#     """
#     try:
#         raw_response = askcos_adapter.askcos_tree_search_raw(input)
#         synth_graph, paths = askcos_tree2synth_paths_with_graph(raw_response, USE_RETRO_RXN_RENDERING=False)
#         return {"graph": SynthGraph(synthesis_graph=synth_graph, target_molecule_node_id="BLANK"), "routes": paths}
#     except Exception as e:
#         raise HTTPException(status_code=420, detail=str(e))

@askcos_wrapper_router.post(
    "/tree-search/parse-paths",
    summary="Get ASKCOS tree search results and parses the routes in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token')",
)
def get_tree_search_parse_paths(input: TreeSearchInput) -> Any:
    """
    Get ASKCOS tree search results and parses the paths in the AICP format, uses direct input from ASKCOS endpoint ('/tree-search/controller/call-sync-without-token').

    Args:
        input (dict): Input object containing the query reaction

    Returns:
        dict: Response object containing the parsed ASKCOS paths
    """
    try:
        raw_response = askcos_adapter.askcos_tree_search_raw_v2(input)
        uds = raw_response["result"]["uds"]
        synth_graph, paths = uds_tree2synth_paths_with_graph(uds, USE_RETRO_RXN_RENDERING=False)
        return {"graph": SynthGraph(synthesis_graph=synth_graph, target_molecule_node_id="BLANK"), "routes": paths}
    except Exception as e:
        raise HTTPException(status_code=420, detail=str(e))


class ConvertToAicpRequest(BaseModel):
    convert_from: Literal["askcos"] = Field(default="askcos", description="The format to convert from", examples=["askcos"])
    source_data: dict = Field(
        ...,
        title="Graph Data",
        description="The graph data to be converted",
        examples=[{"nodes": [], "edges": []}],
    )


# @askcos_wrapper_router.post(
#     "/parse-askcos-tree",
#     summary="Parse ASKCOS tree search results into a SynthGraph",
# )
# def parse_askcos_tree(request: ConvertToAicpRequest) -> Dict:
#     """
#     Parse ASKCOS tree search results into a SynthGraph.

#     Args:
#         input (TreeSearchInput): Input object containing the query reaction

#     Returns:
#         SynthGraph: Parsed synthesis graph
#     """
#     source_data = request.source_data
#     conversion_source = request.convert_from

#     synth_graph, paths = askcos_tree2synth_paths_with_graph(
#         TreeSearchResponse(**source_data), USE_RETRO_RXN_RENDERING=False)

#     # Flatten node and edge data into list of objects
#     nodes = [
#         {**{k: v for k, v in attrs.items() if k != "node_id"}, "node_id": n}
#         for n, attrs in synth_graph.nodes(data=True)
#     ]

#     edges = [
#         dict(source=u, target=v, **attrs)
#         for u, v, attrs in synth_graph.edges(data=True)
#     ]

#     # Loop through routes
#     final_routes = []
#     for route in paths:
#         final_routes.append(
#             {
#                 "aggregated_yield": 0.0,
#                 "predicted": True,
#                 "route_index": route["path_index"],
#                 "route_status": "Predicted Synthesis Route",
#                 "method": "ASKCOS v2",
#                 "route_node_labels": route["nodes"]
#             }
#         )

#     # Return converted graph
#     return {
#         "predicted_synth_graph": {
#             "nodes": nodes,
#             "edges": edges
#         },
#         "routes": final_routes
#     }


# @askcos_wrapper_router.post(
#     "/askcos-condition-predict",
#     summary="Predict reaction conditions using ASKCOS context prediction",
# )
# def askcos_condition_predict(request: AskcosConditionPredictRequest) -> AskcosConditionPredictResponse:
#     """
#     Predicts suitable reaction conditions for a given reaction SMILES using the ASKCOS context prediction model.

#     This endpoint performs the following steps:
#     - Validates and parses the input reaction SMILES.
#     - Extracts reactant and product SMILES.
#     - Calls the ASKCOS adapter to initiate context prediction via an async POST/polling flow.
#     - Returns a structured response containing recommended reaction conditions.

#     Returns:
#         AskcosConditionPredictResponse: The original input rxsmiles and a list of condition options.

#     Raises:
#         HTTPException: If input is malformed, parsing fails, or the ASKCOS service raises an error.
#     """
#     rxsmiles = request.rxsmiles
#     num_results = request.num_results

#     # Basic format validation for reaction SMILES
#     if ">" not in rxsmiles:
#         raise HTTPException(status_code=420, detail="Input must be reaction smiles.")
    
#     # Ensure the reaction SMILES can be parsed by RDKit
#     if not is_reaction_smiles_parseable(rxsmiles):
#         raise HTTPException(status_code=420, detail="Reaction SMILES not parseable by RDKit.")

#     try:
#         # Split the reaction into reactants and products
#         parsed_reaction = decomposition_utils.parse_reaction_smiles(rxsmiles)
#         reactants = ".".join(parsed_reaction.reactants)
#         products = ".".join(parsed_reaction.products)
#     except Exception:
#         raise HTTPException(status_code=420, detail="Failed to parse reaction SMILES.")

#     try:
#         # Call the ASKCOS context prediction API via the adapter
#         response = askcos_adapter.askcos_context_prediction(reactants, products, num_results)
        
#         # Return a structured response with condition options
#         return AskcosConditionPredictResponse(
#             rxsmiles=rxsmiles,
#             condition_options=response.output
#         )
#     except Exception as e:
#         logger.error(f"Error during ASKCOS context prediction: {str(e)}")
#         # Return server error if ASKCOS call fails
#         raise HTTPException(status_code=500, detail=f"Error from ASKCOS context prediction: {str(e)}")    
    