from abc import ABC, abstractmethod

from syngps.models.models import Availability


class AbstractInventoryAdapter(ABC):
    """
    Abstract class for inventory adapter.
    """

    def __init__(self):
        """
        Inits abstract class for inventory adapter.
        """
        super().__init__()

    @abstractmethod
    def in_inventory_by_inchikey(self, inchikey: str) -> bool:
        """
        Method to check if a substance is in inventory by its inchikey.

        Args:
            inchikey (str): Inchikey of the substance.

        Returns:
            bool: True if substance is in inventory, False otherwise.
        """
        pass

    @abstractmethod
    def inchikey_inventory_status(self, inchikeys: list) -> dict:
        """
        Method to get the inventory status of the given inchikeys.

        Args:
            inchikeys (List[str]): List of inchikeys to check inventory status for.

        Returns:
            Dict[str, InventoryStatus]: A dictionary mapping each inchikey to its inventory status.
        """
        pass

    @abstractmethod
    def inventory_status_details(self, inchikey: str) -> Availability:
        """
        Method to get the inventory status for a given inchikey.

        Args:
            inchikey (str): The inchikey to check inventory status for.

        Returns:
            Availability: The availability information for the given inchikey.
        """
        pass

    @abstractmethod
    def adapter_status(self) -> dict:
        """
        Method to get the status of the inventory adapter.

        Returns:
            dict: A dictionary with the connection status to the data sources
        """
        pass
