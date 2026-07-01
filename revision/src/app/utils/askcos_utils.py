from rdkit.Chem import rdChemReactions
from app.models import PredictedRouteDetails, TreeSearchInput


def is_reaction_smiles_parseable(smiles: str) -> bool:
    """
    Checks if a reaction SMILES string is parseable by RDKit.

    Args:
        smiles (str): The reaction SMILES string to evaluate.

    Returns:
        bool: True if the SMILES is a parseable reaction, False otherwise.
    """
    try:
        smiles = smiles.split("|")[0].strip()
        reaction = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
        return reaction is not None
    except Exception:
        return False


def prediction_request2askcos_tree_search_input(target_molecule_smiles: str, reaction_steps: int, input: PredictedRouteDetails) -> TreeSearchInput:
    """
    Converts a PredictedRouteDetails object to an ASKCOS TreeSearchInput object.
    """
    askcos_input = TreeSearchInput(smiles=target_molecule_smiles)
    askcos_input.build_tree_options.expansion_time = 10
    askcos_input.build_tree_options.max_depth = reaction_steps
    askcos_input.enumerate_paths_options.paths_only = True
    askcos_input.enumerate_paths_options.max_paths = input.max_routes
    return askcos_input
