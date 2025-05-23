from abc import ABC, abstractmethod
from typing import List


# Abstract Yield Predicter Class
class AbstractYieldPredictor(ABC):
    """
    Abstract class for yield predicter.
    """

    def __init__(self):
        """
        Inits abstract class for yield predicter.
        """
        super().__init__()

    @abstractmethod
    def predict_reaction_yield(self, rxsmiles: str) -> float:
        """
        Method to predict yield for a given reaction smiles.

        Args:
            rxsmiles (str): Reaction SMILES.

        Returns:
            float: Predicted yield.
        """
        pass

    @abstractmethod
    def predict_multiple_reaction_yield(self, rxsmiles: List[str]) -> List[float]:
        """
        Method to predict yield for a list of reaction smiles.

        Args:
            rxsmiles (List[str]): List of reaction SMILES.

        Returns:
            List[float]: List of predicted yields.
        """
        pass
