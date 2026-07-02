import time
from typing import Any, Optional

import requests
from app.logging_config import logger
from app.models import (
    AskcosAtommapResponse,
    AtommapRequest,
    AtommapResponse,
    ModuleStatusesModel,
    TreeSearchInput,
    TreeSearchResponse,
    ContextPredictionRequest,
    ContextPredictionResponse,
)


class AskcosApiAdapter:
    """
    Adapter class for interacting with the ASKCOS API.

    This class provides methods to send HTTP requests to the ASKCOS API's endpoints,
    handle responses, and parse them into defined models.
    """

    def __init__(self, askcos_base_url: Optional[str] = None, askcos_timeout: int = 120):
        # Initialize the adapter with the ASKCOS base URL. Raises an error if the URL is not provided.
        if not askcos_base_url:
            raise ValueError("askcos_base_url is required (e.g., https://askcos.mit.edu/api).")
        self.askcos_base_url = askcos_base_url
        self.askcos_timeout = askcos_timeout

    def _askcos_post_request(self, endpoint: str, data: dict) -> dict:
        """
        Sends a POST request to a specified ASKCOS API endpoint with the given data payload.

        Args:
            endpoint (str): API endpoint to send the request to.
            data (dict): Payload to send in the body of the POST request.

        Returns:
            dict: JSON response from the ASKCOS API as a dictionary.

        Raises:
            ValueError: If a network error occurs or if the request fails.
        """
        try:
            url = self.askcos_base_url + endpoint  # Construct the full URL
            response = requests.post(url, json=data, timeout=self.askcos_timeout)  # Send POST request with JSON data
            response.raise_for_status()  # Raise an error if the request failed
            return response.json()  # Parse and return the JSON response
        except requests.exceptions.RequestException as e:
            # Log and re-raise the exception as a ValueError with an error message
            logger.error(f"Error connecting to ASKCOS: {e}")
            raise ValueError(f"Error connecting to ASKCOS: {e}")

    def _askcos_get_request(self, endpoint: str, params: Optional[dict] = None) -> dict:
        """
        Sends a GET request to a specified ASKCOS API endpoint with optional query parameters.

        Args:
            endpoint (str): API endpoint to send the request to.
            params (Optional[dict]): Query parameters to include in the request URL.

        Returns:
            dict: JSON response from the ASKCOS API as a dictionary.

        Raises:
            ValueError: If a network error occurs or if the request fails.
        """
        try:
            url = self.askcos_base_url + endpoint  # Construct the full URL
            response = requests.get(url, params=params, timeout=self.askcos_timeout)  # Send GET request with optional query parameters
            response.raise_for_status()  # Raise an error if the request failed
            return response.json()  # Parse and return the JSON response
        except requests.exceptions.RequestException as e:
            # Log and re-raise the exception as a ValueError with an error message
            logger.error(f"Error connecting to ASKCOS: {e}")
            raise ValueError(f"Error connecting to ASKCOS: {e}")

    def get_askcos_status(self) -> ModuleStatusesModel:
        """
        Retrieves the status of ASKCOS backend modules.

        This method calls a specific ASKCOS endpoint to get the status of backend modules,
        checks if all modules are ready, and maps the response into a ModuleStatusesModel object.

        Returns:
            ModuleStatusesModel: The response data mapped into a ModuleStatusesModel instance.
        """
        response = self._askcos_get_request("/admin/get-backend-status")  # Call API to get backend status
        # Check if all modules are 'ready' and store the result in 'all_healthy' key
        all_statuses = [module["ready"] for module in response["modules"]]
        response["all_healthy"] = all(all_statuses)  # 'all_healthy' is True if all modules are ready
        # Map the response dictionary to a ModuleStatusesModel object and return it
        return ModuleStatusesModel(**response)

    def askcos_atommap_rxnmapper(self, request: AtommapRequest) -> AtommapResponse:
        """
        Calls the ASKCOS Atommap endpoint to map atom indices in a reaction SMILES string.

        This method sends a POST request to the ASKCOS Atommap endpoint with the input reaction SMILES string,

        Args:
            request (AtommapRequest): An AtommapRequest object containing the input reaction SMILES string.

        Returns:
            AtommapResponse: An AtommapResponse object containing the mapped reaction SMILES string and confidence.
        """
        try:
            askcos_request = {"smiles": [request.smiles]}

            logger.debug("Making call to ASKCOS '/atom-map/rxnmapper/call-sync' endpoint. Input smiles: " + request.smiles)
            response = AskcosAtommapResponse(**self._askcos_post_request("/atom-map/rxnmapper/call-sync", askcos_request))
            logger.debug("ASKCOS Response code: " + str(response.status_code))

            if response.status_code != 200:
                logger.error(f"Error from ASKCOS Atommap, status {response.status_code}: " + response.message)
                raise Exception("ASKCOS Atommap endpoint returned non 200 status")

            return AtommapResponse(original_rxsmiles=request.smiles, mapped_rxsmiles=response.result[0].mapped_rxn, confidence=response.result[0].confidence)
        except Exception as e:
            logger.error("Error from ASKCOS RXNMapper Atommap: ", exc_info=e)
            raise Exception("Error from ASKCOS RXNMapper Atommap")

    def askcos_tree_search_raw(self, input: TreeSearchInput) -> TreeSearchResponse:
        """
        Calls the ASKCOS Tree Search endpoint to search for reactions based on a query reaction SMILES string. This function takes the
        exact input that ASKCOS offers and will return the raw ASKCOS response.

        Args:
            input (TreeSearchInput): A TreeSearchInput object containing the input reaction SMILES string.

        Returns:
            TreeSearchResponse: A TreeSearchResponse object containing the raw ASKCOS response.
        """
        API_ENDPOINT = "/tree-search/controller/call-sync-without-token"
        try:
            logger.debug(f"Making call to ASKCOS '{API_ENDPOINT}' endpoint. Input:")
            logger.debug(input.model_dump())
            response = self._askcos_post_request(API_ENDPOINT, input.model_dump())
            logger.debug("ASKCOS Response code: " + str(response.get("status_code", "No status code in response")))
            return TreeSearchResponse(**response)
        except requests.exceptions.Timeout:
            logger.error("Timeout error from ASKCOS Tree Search")
            raise Exception("Timeout error from ASKCOS Tree Search")
        except Exception as e:
            logger.error("Error from ASKCOS Tree Search: ", exc_info=e)
            raise Exception("Error from ASKCOS Tree Search")
        
    def askcos_tree_search_raw_v2(self, input: TreeSearchInput) -> Any:
        """
        Calls the ASKCOS Tree Search endpoint to search for reactions based on a query reaction SMILES string. This function takes the
        exact input that ASKCOS offers and will return the raw ASKCOS response, but mapped to the TreeSearchResponseV2 model which has some additional parsing and structuring compared to the original TreeSearchResponse model.

        Args:
            input (TreeSearchInput): A TreeSearchInput object containing the input reaction SMILES string.

        Returns:
            Any: A object containing the raw ASKCOS response mapped to the new model.
        """
        API_ENDPOINT = "/tree-search/controller/call-sync-without-token"
        try:
            logger.debug(f"Making call to ASKCOS '{API_ENDPOINT}' endpoint. Input:")
            logger.debug(input.model_dump())
            response = self._askcos_post_request(API_ENDPOINT, input.model_dump())
            logger.debug("ASKCOS Response code: " + str(response.get("status_code", "No status code in response")))
            return response
        except requests.exceptions.Timeout:
            logger.error("Timeout error from ASKCOS Tree Search")
            raise Exception("Timeout error from ASKCOS Tree Search")
        except Exception as e:
            logger.error("Error from ASKCOS Tree Search: ", exc_info=e)
            raise Exception("Error from ASKCOS Tree Search")
        
    def askcos_context_prediction(self, reactants: str, products: str, num_results: int = 5) -> ContextPredictionResponse:
        CONTEXT_API_ENDPOINT = "/legacy/context"
        ASYNC_TASK_ENDPOINT = "/legacy/celery/task/"

        # Step 1: POST to context endpoint to start task
        context_body = ContextPredictionRequest(
            reactants=reactants,
            products=products,
            num_results=num_results,
            with_smiles=True,
            return_scores=True
        ).model_dump()

        # Send initial request and get task_id
        post_response: dict[str, Any] = self._askcos_post_request(CONTEXT_API_ENDPOINT, context_body)
        task_id = post_response.get("task_id")
        if not task_id:
            raise ValueError("No task_id returned from context prediction request.")

        # Step 2: Poll async task status until SUCCESS or timeout
        timeout = 30  # seconds
        poll_interval = 1  # second
        start_time = time.time()

        while True:
            # Check for timeout
            if time.time() - start_time > timeout:
                raise TimeoutError(f"Task {task_id} did not complete within {timeout} seconds.")

            # Poll task status
            task_status: dict[str, Any] = self._askcos_get_request(ASYNC_TASK_ENDPOINT + task_id)

            state = task_status.get("state")
            complete = task_status.get("complete")
            failed = task_status.get("failed")

            logger.debug(f"Polling task {task_id}: state={state}, complete={complete}, failed={failed}")
            logger.debug(task_status)

            # Check for desired condition
            if state == "SUCCESS" and complete and not failed:
                return ContextPredictionResponse(**task_status)

            time.sleep(poll_interval)

