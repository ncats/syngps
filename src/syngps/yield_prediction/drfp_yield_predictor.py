import os
import pickle
from typing import List

import numpy as np
from drfp import DrfpEncoder
from syngps.errors import ReactionYieldPredictionError
from syngps.yield_prediction import AbstractYieldPredictor
import logging

logger = logging.getLogger(__name__)


# Drfp Yield Predicter Class
class DrfpYieldPredictor(AbstractYieldPredictor):
    """
    Yield predictor class using DRFP method.

    DRFP GitHub: https://github.com/reymond-group/drfp
    DRFP Paper: https://pubs.rsc.org/en/content/articlehtml/2022/dd/d1dd00006c
    """

    def __init__(self):
        """
        Inits predictor class. Loads DRFP model.
        """
        try:
            # Determine the directory of this script
            dir_path = os.path.dirname(os.path.realpath(__file__))

            # Construct the path to the model file relative to the directory of this script
            self.model_path = os.path.join(dir_path, "models", "buchwald_hartwig_model.pkl")

            # Load the model from the constructed path
            with open(self.model_path, "rb") as file:
                self.model = pickle.load(file)
        except Exception:
            logger.error(f"Failed to load model: {self.model_path}", exc_info=True)
            raise ValueError(f"Failed to load model: {self.model_path}")

    def predict_reaction_yield(self, rxsmiles: str) -> float:
        """
        Method to predict yield for a given reaction smiles.

        Args:
            rxsmiles (str): Reaction SMILES.

        Returns:
            float: Predicted yield.
        """
        try:
            X, mapping = DrfpEncoder.encode(
                rxsmiles,
                n_folded_length=2048,
                radius=3,
                rings=True,
                mapping=True,
            )

            X = np.asarray(
                X,
                dtype=np.float32,
            )
            # Following comment used to be a parameter, no longer accepted
            # ntree_limit=self.model.best_ntree_limit
            yield_prediction = self.model.predict(X)
        except Exception:
            logger.error(f"Failed to predict yield for reaction: {rxsmiles}", exc_info=True)
            raise ReactionYieldPredictionError(rxsmiles, "Failed drfp prediction")
        return yield_prediction[0]

    def predict_multiple_reaction_yield(self, rxsmiles: List[str]) -> List[float]:
        """
        Method to predict yield for a list of reaction smiles.

        Args:
            rxsmiles (List[str]): List of reaction SMILES.

        Returns:
            List[float]: List of predicted yields.
        """
        try:
            yields = [self.predict_reaction_yield(rxsmiles) for rxsmiles in rxsmiles]
        except ReactionYieldPredictionError as e:
            logger.error(f"Failed to predict yield for reactions: {rxsmiles}", exc_info=True)
            raise e
        return yields
