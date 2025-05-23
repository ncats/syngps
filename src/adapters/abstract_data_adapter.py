from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Set

from models import SynthGraph, SynthGraphSearch


# Abstract Data Adapter Class
class AbstractDataAdapter(ABC):
    """
    Abstract class for data adapter.
    """

    def __init__(self):
        """
        Inits abstract class for data adapter.
        """
        super().__init__()

    @abstractmethod
    def verify_connection(self) -> bool:
        """
        Method to verify connection to the knowledge base.

        Returns:
            bool: True if connection is successful, False otherwise.
        """
        pass

    @abstractmethod
    def close(self):
        """
        Method to close connection to the knowledge base.
        """
        pass

    @abstractmethod
    def get_data_counts(self) -> Dict[str, int]:
        pass

    @abstractmethod
    def get_reaction_by_rxid(self, rxid: str) -> Dict[str, Any]:
        """
        Method to get reaction for a given reaction ID.

        Args:
            rxid (str): Reaction ID.

        Returns:
            Optional[Dict[str, Any]]: Reaction details if found, None otherwise.
        """
        pass

    @abstractmethod
    def get_rxsmiles_extended(self, rxid: str) -> str:
        """
        Method to get extended rxsmiles for a given reaction ID.

        Args:
            rxid (str): Reaction ID.

        Returns:
            str: Extended rxsmiles if found, None otherwise.
        """
        pass

    @abstractmethod
    def get_multiple_rxsmiles_extended(self, rxids: List[str]) -> Dict[str, Optional[str]]:
        """
        Method to get multiple extended rxsmiles for a list of reaction IDs.

        Args:
            rxids (List[str]): List of reaction IDs.

        Returns:
            Dict[str, Optional[str]]: Dictionary of rxid to extended rxsmiles.
        """
        pass

    @abstractmethod
    def find_target_molecule_node(self, target_molecule: str) -> str:
        """
        Method to find the target molecule node.

        Args:
            target_molecule (str): Target molecule.

        Returns:
            str: Node ID of the target molecule.
        """
        pass

    @abstractmethod
    def fetch_ambiguous_reactions(self) -> Set[str]:
        """
        Method to find ambiguous reactions in the knowledge base.

        Returns:
            Set[str]: Set of ambiguous reactions.
        """
        pass

    @abstractmethod
    def fetch_multiproduct_reactions(self) -> Set[str]:
        """
        Method to find multiproduct reactions in the knowledge base.

        Returns:
            Set[str]: Set of multiproduct reactions.
        """
        pass

    @abstractmethod
    def fetch_synthesis_graph(self, search_params: SynthGraphSearch) -> SynthGraph:
        """
        Method to fetch the synthesis graph based on search parameters.

        Args:
            search_params (SynthGraphSearch): Search parameters.

        Returns:
            SynthGraph: Synthesis graph.
        """
        pass

    @abstractmethod
    def get_substance_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Method to get substance for a given InChIKey.

        Args:
            inchikey (str): InChIKey.

        Returns:
            Optional[Dict[str, Any]]: Substance details if found, None otherwise.
        """
        pass
