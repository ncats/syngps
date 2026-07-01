import json
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

from fastapi import APIRouter, Body, HTTPException
from networkx import DiGraph
from app.adapters import KBAdapter, AskcosInventoryAdapter, CustomStockInventoryAdapter
from syngps import (
    AggregateYieldsInput,
    AggregateYieldsOutput,
    AicpFunctions,
    InchikeyNotFoundError,
    ParseSynthGraphInput,
    ParseSynthGraphOutput,
    ReactionNode,
    SGPInput,
    SynthGraph,
    SynthGraphSearch,
    SynthRoute,
    SynthRoutePreDefined,
    TopNYieldSearch,
    TopNYieldSynthRoutesResult,
)
from syngps.errors import (
    SubstanceNotFoundInSynthGraphError,
    SynthGraphParsingException,
)
from syngps.utils import (
    smiles2inchikey,
    json_routes_to_graph_subgraphs,
    smiles2inchikey,
)
from syngps.utils.graph_utils import merge_synth_graphs
from app.app_config import APP_CONFIG
from app.logging_config import logger
from app.models import (
    NoPathsFoundInAskcosResponse,
    NoResultFoundInAskcosResponse,
    SynthesisRoutesRequest,
    SynthesisRoutesResponse,
)
from app.routers.asckos_wrappers import (
    get_askcos_status,
    get_tree_search_raw,
)
from app.utils import uds_tree2synth_paths_with_graph
from app.utils import prediction_request2askcos_tree_search_input

synthplanning_router = APIRouter(prefix="/prediction", tags=["Synthplanning"])

adapter = KBAdapter(
    mongo_uri=APP_CONFIG.mongo_connection_uri,
    mongo_db_name=APP_CONFIG.mongo_db_name,
    mongo_user=APP_CONFIG.mongo_db_user,
    mongo_password=APP_CONFIG.mongo_db_pass,
    memgraph_uri=APP_CONFIG.memgraph_connection_uri,
    memgraph_user=APP_CONFIG.memgraph_user,
    memgraph_password=APP_CONFIG.memgraph_pass,
    memgraph_conn_encrypted=APP_CONFIG.memgraph_conn_encrypted,
    smiles_encryption_key=APP_CONFIG.smiles_encryption_key,
    graph_backend="memgraph",
)

# Set up inventory adapters, require all to be defined
if not all([APP_CONFIG.askcos_inv_file_path, APP_CONFIG.custom_stock_inv_file_path]):
    raise ValueError("All inventory file paths (ASKCOS_INV_FILE_PATH, CUSTOM_STOCK_INV_FILE_PATH) must be defined in the configuration.")

with ThreadPoolExecutor(max_workers=4) as executor:
    askcos_future = executor.submit(AskcosInventoryAdapter, askcos_inv_file_path=APP_CONFIG.askcos_inv_file_path)
    custom_stock_future = executor.submit(CustomStockInventoryAdapter, stock_inv_file_path=APP_CONFIG.custom_stock_inv_file_path)

askcos_inv_adapter = askcos_future.result()
custom_stock_adapter = custom_stock_future.result()

aicp_with_askcos_synthplanning = AicpFunctions(data_adapter=adapter, inventory_adapter=askcos_inv_adapter)
aicp_with_custom_stock_synthplanning = AicpFunctions(data_adapter=adapter, inventory_adapter=custom_stock_adapter)


SERVICE_NAME = "synthplanning"

######################
# Custom data models #
######################


######################
#     Endpoints      #
######################
@synthplanning_router.get("/status", summary="Get the status of the synthplanning service and relevant data sources")
def synthplanning_status() -> dict:
    """
    Validate the connection to the data sources.

    Returns:
        dict: A dictionary with the connection status to the data sources
    """
    return {
        "service": SERVICE_NAME,
        "mongo_connected": adapter.verify_mongo_connection(),
        "graphdb_connected": adapter.verify_graphdb_connection(),
        "askcos_connected": get_askcos_status().all_healthy,
        "askcos_inventory_adapter": askcos_inv_adapter.adapter_status(),
        "data": adapter.get_data_counts(),
    }


@synthplanning_router.post(
    "/fetch_synthesis_graph",
    summary="Fetch the synthesis graph based on the search parameters",
    response_model_exclude_none=True,
    response_model_exclude_unset=True,
)
async def fetch_synthesis_graph(search_params: SynthGraphSearch) -> SynthGraph:
    """
    Fetch the synthesis graph for a given target molecule InChIKey using the specified search parameters.

    **Arguments**:
    - `target_molecule_inchikey` (str): InChIKey of the target molecule.
      Example: `"YTIQRXMAVJLXHH-UHFFFAOYSA-N"`
    - `reaction_steps` (int): Number of reaction steps to search away from the target molecule.
      Default: `2`
    - `query_type` (str): Type of query used to extract graph data.
      Default: `"shortest_path"`, which uses shortest path logic for faster queries.
    - `leaves_as_sm` (bool): Whether to treat leaf nodes in the graph as starting materials.
      Default: `True`
    - `include_availability_info` (bool): Whether to include availability information in the response.
      Default: `False`
    - `annotate_reactions` (bool): Whether to annotate reactions in the graph with Hazelnut.
      Default: `False`
    - `inventory_source` (str): Source of inventory to check for availability.
      Default: `"askcos"`, options are `"askcos"`, `"stock"` (for custom stock inventory).

    **Returns**:
    - `SynthGraph`: A dictionary containing the synthesis graph data.
    """

    try:
        # Select the appropriate AICP synthplanning function based on the inventory source specified in the search parameters
        if search_params.inventory_source == "askcos":
            aicp_synthplanning = aicp_with_askcos_synthplanning
        elif search_params.inventory_source == "stock":
            aicp_synthplanning = aicp_with_custom_stock_synthplanning
        else:
            raise HTTPException(status_code=400, detail=f"Invalid inventory source: {search_params.inventory_source}. Valid options are 'askcos', 'stock'.")
        # Fetch the synthesis graph based on the search parameters
        graph = aicp_synthplanning.fetch_synthesis_graph(search_params)
        return graph
    except HTTPException as e:
        raise e
    except SubstanceNotFoundInSynthGraphError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Error fetching synthesis graph: {e}")
        raise HTTPException(status_code=500, detail="Error fetching synthesis graph")


@synthplanning_router.post(
    "/top_yield_routes",
    summary="Top N synthesis routes based on overall yield, using only evidence based reactions from the synthesis graph",
    response_model_exclude_none=True,
    response_model_exclude_unset=True,
)
def top_yield_routes(search_params: TopNYieldSearch) -> TopNYieldSynthRoutesResult:
    """
    Retrieves the top N synthesis routes based on the overall yield.
    Fetches the synthesis graph for a given target molecule InChIKey using the specified search parameters.

    **Arguments**:
    - `target_molecule_inchikey` (str): InChIKey of the target molecule.
      Example: `"YTIQRXMAVJLXHH-UHFFFAOYSA-N"`
    - `reaction_steps` (int): Number of reaction steps to search away from the target molecule.
      Default: `2`
    - `query_type` (str): Type of query used to extract graph data.
      Default: `"shortest_path"`, which uses shortest path logic for faster queries.
    - `leaves_as_sm` (bool): Whether to treat leaf nodes in the graph as starting materials.
      Default: `True`
    - `graph_backend` (str): Graph database backend to use, such as `"memgraph"` or `"neo4j"`.
      Default: `"memgraph"`
    - `top_n_routes` (int): Number of top synthesis routes to return.
      Default: `3`
    - `include_route_candidates` (bool): Whether to include route candidate structures in the result.
      Default: `False`
    - `include_combination_graphs` (bool): Whether to include combination graphs in the output.
      Default: `False`
    - `synthesis_graph_json` (Optional[Dict[str, Any]]): Optionally provide an existing synthesis graph as input instead of querying from the database.
      Default: `None`
    - `inventory_source` (str): Source of inventory to check for availability, options are 'askcos', 'stock'. Default is 'askcos'.


    **Returns**:
    - `TopNYieldSynthRoutesResult` object - contains the top N synthesis routes based on the overall yield
    """

    try:
        # Select the appropriate AICP synthplanning function based on the inventory source specified in the search parameters
        if search_params.inventory_source == "askcos":
            aicp_synthplanning = aicp_with_askcos_synthplanning
        elif search_params.inventory_source == "stock":
            aicp_synthplanning = aicp_with_custom_stock_synthplanning
        else:
            raise HTTPException(status_code=400, detail=f"Invalid inventory source: {search_params.inventory_source}. Valid options are 'askcos', 'stock'.")
        
        result = aicp_synthplanning.find_top_n_yield_synthesis_routes(search_params, top_n=search_params.top_n_routes)

        # Remove certain fields from the response to reduce the payload size
        result.synth_graph.search_params = None

        # Remove combination graphs and route candidates if not requested
        if not search_params.include_combination_graphs:
            result.combination_graphs = None
        if not search_params.include_route_candidates:
            result.route_candidates = None

        # Move availability from the synthesis graph to the result
        if search_params.include_availability_info:
            result.availability = result.synth_graph.availability
            result.synth_graph.availability = None

        return result
    except SynthGraphParsingException as e:
        raise HTTPException(status_code=400, detail=str(e))
    except SubstanceNotFoundInSynthGraphError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        import traceback

        stack_trace = traceback.format_exc()
        logger.error(f"Error fetching top n routes: {e}\nStacktrace:\n{stack_trace}")
        raise HTTPException(status_code=500, detail="Error fetching top n routes")


@synthplanning_router.post(
    "/synthesis_routes",
    summary="Get synthesis routes for the given target molecule. Combines both evidence-based and predicted routes.",
    response_model_exclude_none=True,
    response_model_exclude_unset=True,
)
def synthesis_routes(request: SynthesisRoutesRequest) -> SynthesisRoutesResponse:
    """
    Retrieve synthesis routes for a given target molecule.

    Either a target SMILES (`target_molecule_smiles`) or an InChIKey (`target_molecule_inchikey`) must be provided.
    If both are provided, they must correspond to the same molecule.

    Either `include_evidence_routes` or `include_predicted_routes` must be set to `True`.

    **Arguments**:
    - `target_molecule_inchikey` (str): InChIKey of the target molecule
    - `target_molecule_smiles` (str): SMILES of the target molecule
    - `reaction_steps` (int): Maximum number of steps in the synthesis
    - `include_svgs` (bool): Whether to include SVGs in the response
    - `include_evidence_routes` (bool): Include evidence-based routes if True
    - `evidence_options` (EvidenceBasedRouteDetails): Parameters for evidence-based route generation
    - `include_evidence_synth_graph` (bool): Include evidence-based synthesis graph if True
    - `include_predicted_routes` (bool): Include predicted routes if True
    - `prediction_options` (PredictedRouteDetails): Parameters for predictive route generation
    - `include_predicted_synth_graph` (bool): Include predictive synthesis graph if True
    - `inventory_source` (str): Source of inventory to check for availability, options are 'askcos', 'stock'. Default is 'askcos'.

    **Returns**:
    `SynthesisRoutesResponse`: The response containing synthesis routes
    """
    target_smiles, target_inchikey = _resolve_target_molecule(request, adapter)

    request.target_molecule_inchikey = target_inchikey
    request.target_molecule_smiles = target_smiles

    include_evidence_routes = request.include_evidence_routes
    include_predicted_routes = request.include_predicted_routes

    if not include_evidence_routes and not include_predicted_routes:
        raise HTTPException(
            status_code=400,
            detail="At least one of include_evidence_routes or include_predicted_routes must be True",
        )

    try:
        with ThreadPoolExecutor() as executor:
            futures: Dict[Future, str] = {}

            if include_evidence_routes:
                futures[
                    executor.submit(
                        _get_evidence_routes,
                        TopNYieldSearch(
                            target_molecule_inchikey=target_inchikey,
                            reaction_steps=request.reaction_steps,
                            query_type=request.evidence_options.query_type,
                            leaves_as_sm=True,
                            inventory_source=request.inventory_source,
                        ),
                        top_n=request.evidence_options.top_n_routes,
                    )
                ] = "evidence"
            else:
                evidence_based_success = True

            if include_predicted_routes:
                futures[
                    executor.submit(
                        _get_predictive_routes,
                        target_smiles,
                        request.reaction_steps,
                        request.prediction_options,
                    )
                ] = "predictive"
            else:
                predictive_success = True

            top_n_routes_result = None
            predictive_routes = None

            # Wait for both tasks to complete and collect results
            for future in as_completed(futures):
                result_type = futures[future]
                if result_type == "evidence":
                    evidence_based_success, top_n_routes_result = future.result()
                elif result_type == "predictive":
                    predictive_success, predictive_graph, predictive_routes = future.result()

        if (include_evidence_routes and  not evidence_based_success) and (include_predicted_routes and not predictive_success):
            raise HTTPException(status_code=500, detail="Failed to retrieve both evidence based and predictive synthesis routes")

        merged_routes = _merge_routes(top_n_routes_result, predictive_routes)

        if not merged_routes:
            raise HTTPException(status_code=404, detail="No synthesis routes found")

        # Prepare graphs for merging
        evidence_graph = (
            top_n_routes_result.synth_graph.synthesis_graph
            if (
                include_evidence_routes and
                evidence_based_success
                and top_n_routes_result is not None
                and hasattr(top_n_routes_result, "synth_graph")
                and top_n_routes_result.synth_graph is not None
            )
            else None
        )
        predictive_graph_obj = predictive_graph if (predictive_success and predictive_graph is not None) else None

        merged_synth_graph = None
        if evidence_graph and predictive_graph_obj:
            merged_synth_graph = merge_synth_graphs([evidence_graph, predictive_graph_obj], target_inchikey)

        # Merge the routes and return
        response = SynthesisRoutesResponse(
            target_molecule_inchikey=target_inchikey,
            target_molecule_smiles=target_smiles,
            reaction_steps=request.reaction_steps,
            evidence_routes_success=evidence_based_success,
            predicted_routes_success=predictive_success,
            routes=merged_routes,
        )

        # Optionally add evidence synthesis graph
        if evidence_based_success and request.include_evidence_synth_graph and top_n_routes_result is not None:
            response.evidence_synth_graph = top_n_routes_result.synth_graph

        # Optionally add predicted synthesis graph
        if predictive_success and request.include_predicted_synth_graph and predictive_graph is not None:
            response.predicted_synth_graph = SynthGraph(target_molecule_node_id=target_inchikey, synthesis_graph=predictive_graph, search_params=None)

        # Optionally add merged synthesis graph
        if merged_synth_graph is not None:
            response.merged_synth_graph = SynthGraph(target_molecule_node_id=target_inchikey, synthesis_graph=merged_synth_graph, search_params=None)

        return response
    except SubstanceNotFoundInSynthGraphError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except HTTPException as e:
        raise e
    except Exception:
        logger.error("Error fetching synthesis routes:", exc_info=True)
        raise HTTPException(status_code=500, detail="Error fetching synthesis routes")


def _get_predictive_routes(target_molecule_smiles: str, reaction_steps: int, prediction_options: Any) -> Tuple[bool, DiGraph, List[Dict[str, Any]]]:
    """
    Retrieves the predictive synthesis routes for the given substances.
    """
    try:
        logger.info("Creating ASKCOS input for predictive synthesis routes")
        askcos_input = prediction_request2askcos_tree_search_input(
            target_molecule_smiles=target_molecule_smiles, reaction_steps=reaction_steps, input=prediction_options
        )
        logger.info(f"ASKCOS Search criteria: {askcos_input.model_dump_json(indent=2)}")

        raw_response = get_tree_search_raw(askcos_input)
        uds = raw_response["result"]["uds"]
        predicted_graph, predicted_routes = uds_tree2synth_paths_with_graph(uds)
        return (True, predicted_graph, predicted_routes)
    except NoResultFoundInAskcosResponse:
        return (False, None, [])  # Error in askcos response
    except NoPathsFoundInAskcosResponse:
        return (True, None, [])  # No askcos routes found succesfully
    except Exception:
        logger.error("Error fetching predictive synthesis routes:", exc_info=True)
        return (False, None, [])


def _get_evidence_routes(search_params: TopNYieldSearch, top_n: int) -> Tuple[bool, TopNYieldSynthRoutesResult | None]:
    """
    Retrieves the evidence based synthesis routes for the given substances.
    """
    if search_params.inventory_source == "askcos":
        aicp_synthplanning = aicp_with_askcos_synthplanning
    elif search_params.inventory_source == "stock":
        aicp_synthplanning = aicp_with_custom_stock_synthplanning
    else:
        raise HTTPException(status_code=400, detail=f"Invalid inventory source: {search_params.inventory_source}. Valid options are 'askcos', 'stock'.")

    try:
        results = aicp_synthplanning.find_top_n_yield_synthesis_routes(search_params, top_n=top_n)
        return (True, results)
    except InchikeyNotFoundError:
        return (False, None)
    except SubstanceNotFoundInSynthGraphError:
        return (False, None)
    except Exception:
        logger.error("Error fetching evidence based synthesis routes:", exc_info=True)
        return (False, None)


def _merge_routes(
    top_n_routes_result: TopNYieldSynthRoutesResult | None, predictive_routes: List[Dict[Any, Any]] | None
) -> List[Union[SynthRoute, SynthRoutePreDefined]]:
    """
    Merges the evidence based and predictive synthesis routes.
    """
    idx = 1
    routes: List[Union[SynthRoute, SynthRoutePreDefined]] = []

    # Convert evidence based routes to SynthRoute objects
    if top_n_routes_result is not None and top_n_routes_result.routes is not None:
        if top_n_routes_result.routes is not None and len(top_n_routes_result.routes) > 0:
            for route in top_n_routes_result.routes:
                route.route_index = idx
                route.predicted = False
                routes.append(route)
                idx += 1

    if predictive_routes is not None and len(predictive_routes) > 0:
        # Convert predictive routes to SynthRoutePreDefined objects
        for proute in predictive_routes:
            routes.append(
                SynthRoutePreDefined(
                    route_index=idx,
                    route_node_labels=proute["nodes"],
                    aggregated_yield=None,
                    predicted=True,
                    source="ASKCOS v2",
                )
            )
            idx += 1

    return routes


def _resolve_target_molecule(request, adapter) -> Tuple[str, str]:
    target_smiles = request.target_molecule_smiles
    target_inchikey = request.target_molecule_inchikey

    if target_smiles == "":
        target_smiles = None

    if target_inchikey == "":
        target_inchikey = None

    if not target_smiles and not target_inchikey:
        raise HTTPException(status_code=422, detail="Either target_molecule_smiles or target_molecule_inchikey must be provided")

    if target_smiles and target_inchikey:  # Verify they match
        try:
            computed_inchikey = smiles2inchikey(target_smiles)
        except Exception:
            raise HTTPException(status_code=400, detail=f"Failed to convert SMILES to InChIKey: {target_smiles}")
        if computed_inchikey != target_inchikey:
            raise HTTPException(status_code=400, detail="The provided SMILES and InChIKey do not correspond to the same molecule")
    elif target_smiles:
        try:
            target_inchikey = smiles2inchikey(target_smiles)
        except Exception:
            raise HTTPException(status_code=400, detail=f"Failed to convert SMILES to InChIKey: {target_smiles}")
    elif target_inchikey and len(target_inchikey) > 0:
        try:
            target_smiles = adapter.get_substance_by_inchikey(target_inchikey)["canonical_smiles"]
        except InchikeyNotFoundError:
            raise HTTPException(status_code=400, detail="Could not determine the canonical SMILES for the target molecule, please provide the SMILES directly")
    else:
        raise HTTPException(status_code=422, detail="Either target_molecule_smiles or target_molecule_inchikey must be provided")

    if not target_smiles or not target_inchikey:
        raise HTTPException(status_code=500, detail="Error fetching target molecule information")

    return target_smiles, target_inchikey
