from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from app.app_config import APP_CONFIG
from app.logging_config import logger
from app.models import InventoryStatusRequest, InventoryStatusResponse

###########################
# Inventory Status Router
###########################
inventory_router_v2 = APIRouter(prefix="/inventory", tags=["Synthplanning"])

@inventory_router_v2.get("/asi_status", summary="Validate the connection to the data sources")
def validate_data_source() -> JSONResponse:
    """
    Validate the connection to the data sources.

    Returns:
        dict: A dictionary with the connection status to the data sources
    """
    try:
        response = asi_adapter.asi_adapter_status()
        if not response["ready"] or not response["asi_online"]:
            # Return response body with status code 503
            return JSONResponse(content=response, status_code=503)
        return JSONResponse(content=response, status_code=200)
    except Exception:
        logger.exception("Error retrieving ASI connection status", exc_info=True)
        raise HTTPException(status_code=500, detail="Error retrieving ASI connection status")


@inventory_router_v2.post(
    "/inventory_status", summary="Get the inventory status of the given substances", response_model_exclude_none=True, response_model_exclude_unset=True
)
def inventory_status(request: InventoryStatusRequest) -> InventoryStatusResponse:
    """
    Check the inventory status of the given substances against ASI. The inventory status is a dictionary with the connection status to the data sources.

    Available statuses:
     * IN_STOCK ("In Stock - Stereo Match") - Exact Inchikey match in ASI
     * IN_STOCK_NS ("In Stock - Non-Stereo Match") - NS Inchikey match in ASI
     * NOT_AVAILABLE ("Not Available") - Not available in ASI
     * UNKNOWN ("Unknown") - Unknown ASI status

    Returns:
        InventoryStatusResponse: The inventory status of the given substances
    """
    try:
        inventory_statuses = asi_adapter.inchikey_inventory_status(request.inchikeys)
        return InventoryStatusResponse(inchikeys=request.inchikeys, inventory_statuses=inventory_statuses)
    except Exception:
        logger.exception("Error retrieving inventory status", exc_info=True)
        raise HTTPException(status_code=500, detail="Error retrieving inventory status")
