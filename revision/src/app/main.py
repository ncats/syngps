import asyncio
import uuid

from fastapi import FastAPI
from fastapi import __version__ as fastapi_version
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, RedirectResponse
from app.app_config import APP_CONFIG
from app.logging_config import TID_CONTEXT_VAR, logger
from app.routers import (
    askcos_adapter,
    askcos_wrapper_router,
    inventory_router_v2,
    knowledgebase_adapter,
    reaction_utils_router,
    substance_utils_router,
    synthplanning_router,
)

logger.info(f"Starting FastAPI app on port {APP_CONFIG.app_port}...")

app_prefix = APP_CONFIG.app_prefix.strip("/")
app_prefix = f"/{app_prefix}" if app_prefix else ""

# Defien FastAPI app
app = FastAPI(
    description="Synthplanning Service",
    version="0.0.1",
    servers=None,
    docs_url=f"{app_prefix}/docs",
    redoc_url=f"{app_prefix}/redoc",
    openapi_url=f"{app_prefix}/openapi.json",
    swagger_ui_oauth2_redirect_url=f"{app_prefix}/docs/oauth2-redirect",
)

# CORS Middleware
# TODO - Define origins from configuration
origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Add middleware to add transaction ID to request context
@app.middleware("http")
async def add_transaction_id(request, call_next):
    tid = request.headers.get("AICP-TID", None)
    if not tid:
        tid = str(uuid.uuid4())
    TID_CONTEXT_VAR.set(tid)
    response = await call_next(request)
    return response


# Redirect base route to FastAPI docs
@app.get("/", include_in_schema=False)
async def root_redirect():
    return RedirectResponse(url=f"{app_prefix}/docs")


# Define root endpoint
@app.get(f"{app_prefix}/api/v1", tags=["App"])
async def app_root():
    return {"app": "synthplanning", "app_version": APP_CONFIG.app_version, "fast_api_version": fastapi_version}


# Define main app status endpoint for Nagios, returns 500 if not all services are connected
@app.get(f"{app_prefix}/api/v1/status", tags=["App"])
async def app_status():
    """
    Runs all verification checks concurrently and returns the results.
    If any verification check fails, returns a 400 status.
    """
    ASKCOS_STATUS_TIMEOUT = 10  # seconds — ASKCOS can be slow; don't block the whole status check

    async def check_askcos():
        try:
            return await asyncio.wait_for(
                asyncio.to_thread(lambda: askcos_adapter.get_askcos_status().all_healthy),
                timeout=ASKCOS_STATUS_TIMEOUT,
            )
        except (asyncio.TimeoutError, Exception):
            return False

    # Define all verification functions
    tasks = {
        "graphdb_status": asyncio.to_thread(knowledgebase_adapter.verify_graphdb_connection, quick_query=True),
        "mongodb_status": asyncio.to_thread(knowledgebase_adapter.verify_mongo_connection),
        "askcos_connection_status": check_askcos(),
    }

    # Run all tasks concurrently and gather results
    results = await asyncio.gather(*tasks.values(), return_exceptions=True)

    # Map results back to task names
    status_results = {task_name: (result if not isinstance(result, Exception) else False) for task_name, result in zip(tasks.keys(), results)}

    # Determine overall status
    all_healthy = all(status_results.values())

    # Prepare the response JSON
    response = {
        "app": "synthplanning",
        "app_version": APP_CONFIG.app_version,
        "fast_api_version": fastapi_version,
        "status": "healthy" if all_healthy else "unhealthy",
        "services": status_results,
    }

    # Return 500 if any service is unhealthy
    if not all_healthy:
        return JSONResponse(content=response, status_code=500)
    return JSONResponse(content=response, status_code=200)


#####################
# Define V1 routers
#####################

# Include synthplanning router
app.include_router(synthplanning_router, prefix=f"{app_prefix}/api/v1")

# Include inventory router
app.include_router(inventory_router_v2, prefix=f"{app_prefix}/api/v1")
# Include reaction utils router
app.include_router(reaction_utils_router, prefix=f"{app_prefix}/api/v1")

# Include substance utils router
app.include_router(substance_utils_router, prefix=f"{app_prefix}/api/v1")

# Include ASKCOS wrappers router
app.include_router(askcos_wrapper_router, prefix=f"{app_prefix}/api/v1")
