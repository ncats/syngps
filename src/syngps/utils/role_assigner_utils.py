# Author: Gergely Zahoranszky-Kohalmi PhD, Nathan Miller
#
# Email: gergely.zahoranszky-kohalmi@nih.gov, nathan.miller@nih.gov
#
# Organization: NCATS/NIH
#
# Aim: Distinguish substances into reactants and reagents based on atom mapping
#
# If a substance is on the left hand-side of the reaction equation but none of its atoms are found in the product, it's a reagent
# otherwise it's a reactant.
#
# REMARK
#
# This code is adapted from the original code written by Gergely. This code can be found in the 'references' folder at the root of
# the repository. Additionally some code was adapted as part of 'decomposition_utils'.
#
from typing import Any, List

from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms
from syngps.errors import RxsmilesAtomMappingException
from syngps.utils import decomposition_utils
import logging

logger = logging.getLogger(__name__)



def determine_max_reactant_index(rxsmiles: str) -> int:
    # Extract the SMILES from the rxsmiles if it has an extension
    if decomposition_utils.has_extension(rxsmiles):
        rxsmiles = decomposition_utils.extract_smiles(rxsmiles)

    # Split the rxsmiles to extract reactants and count the number of fragments
    reactants = rxsmiles.split(">")[0].strip().split(".")
    return len(reactants) - 1


def determine_max_reactant_index_from_parsed_rxn(parsed_rxn: decomposition_utils.ParsedReaction) -> int:
    return determine_max_reactant_index(parsed_rxn.rxsmiles)


def determine_min_product_index(rxsmiles: str) -> int:
    # Extract the SMILES from the rxsmiles if it has an extension
    if decomposition_utils.has_extension(rxsmiles):
        rxsmiles = decomposition_utils.extract_smiles(rxsmiles)

    # Split the rxsmiles into reactants, reagents, and products
    parts = rxsmiles.split(">")
    reactants_count = len(parts[0].strip().split("."))
    reagents_count = len(parts[1].strip().split(".")) if parts[1].strip() else 0

    # Minimum product index is the sum of reactants and reagents
    return reactants_count + reagents_count


def determine_min_product_index_from_parsed_rxn(parsed_rxn: decomposition_utils.ParsedReaction) -> int:
    return determine_min_product_index(parsed_rxn.rxsmiles)


def determine_max_product_index(rxsmiles: str) -> int:
    # Extract the SMILES from the rxsmiles if it has an extension
    if decomposition_utils.has_extension(rxsmiles):
        rxsmiles = decomposition_utils.extract_smiles(rxsmiles)

    # Split the rxsmiles into reactants, reagents, and products
    parts = rxsmiles.split(">")
    reactants_count = len(parts[0].strip().split("."))
    reagents_count = len(parts[1].strip().split(".")) if parts[1].strip() else 0
    products_count = len(parts[2].strip().split("."))

    # Maximum product index includes all reactants, reagents, and products
    return reactants_count + reagents_count + products_count - 1


def determine_max_product_index_from_parsed_rxn(parsed_rxn: decomposition_utils.ParsedReaction) -> int:
    return determine_max_product_index(parsed_rxn.rxsmiles)


def are_fragment_indices_coherent(fragment_block, max_reactant_index, min_product_index, max_product_index):
    first = True
    first_index_type = None
    index_type = None

    for fragment in fragment_block:

        if int(fragment) <= int(max_reactant_index) and int(fragment) >= int(0):
            index_type = "reactant"
        elif int(fragment) >= int(min_product_index) and int(fragment) <= int(max_product_index):
            index_type = "product"
        elif int(fragment) > int(max_reactant_index) and int(fragment) < int(min_product_index):
            index_type = "reagent"
        else:
            error_message = "[ERROR] Not all original fragment indices are unique in the extension. Please check RXSMILES."
            logger.error(error_message)
            raise decomposition_utils.FragmentGroupError(error_message)

        if first:
            first_index_type = index_type
            first = False
        else:
            if index_type != first_index_type:
                return False

    return True


def validate_f_section_indices(rxsmiles: str) -> bool:
    parsed_rxn = decomposition_utils.parse_reaction_smiles(rxsmiles)

    if not parsed_rxn.fragment_groups:
        raise ValueError(f"[ERROR] No fragment groups found in the parsed RXSMILES: {rxsmiles}")

    # Only these are defined as the reagents are either present of not
    # Always true: min_reactant_index = 0    # zero-indexed
    max_reactant_index = determine_max_reactant_index_from_parsed_rxn(parsed_rxn)  # zero-indexed
    min_product_index = determine_min_product_index_from_parsed_rxn(parsed_rxn)  # zero-indexed
    max_product_index = determine_max_product_index_from_parsed_rxn(parsed_rxn)  # zero-indexed

    for fragment_block in parsed_rxn.fragment_groups:
        if not are_fragment_indices_coherent(fragment_block, max_reactant_index, min_product_index, max_product_index):
            logger.error("[ERROR] Role section border indexing error. Please check the code logic.")
            raise decomposition_utils.FragmentGroupError("[ ERROR] Role section border indexing error. Please check the code logic.")

    # Check uniquness of fragment indices, each must only appear in one block
    unique_fragment_indices: dict[str | int, bool] = {}

    for fragment_block in parsed_rxn.fragment_groups:
        for fragment_index in fragment_block:
            if fragment_index not in unique_fragment_indices.keys():
                unique_fragment_indices[fragment_index] = True
            else:
                logger.error(
                    f"[ERROR] Not all original fragment indices are unique in the extension. Please check RXSMILES. Skipping RXSMILES {rxsmiles}, {parsed_rxn.fragment_groups}."
                )
                raise decomposition_utils.FragmentGroupError(
                    f"[ERROR] Not all original fragment indices are unique in the extension. Please check RXSMILES. Skipping RXSMILES {rxsmiles}, {parsed_rxn.fragment_groups}."
                )

    return True


def rxsmiles_has_atommapping(rxsmiles: str) -> bool:
    """
    Checks if a given reaction SMILES string includes atom mapping.

    Args:
        rxsmiles (str): A reaction SMILES string in the format 'reactants > reagents > products'.

    Returns:
        bool: True if atom mapping is found, False otherwise.
    """
    try:
        # Parse the reaction SMILES into a reaction object
        rxn = AllChem.ReactionFromSmarts(rxsmiles)
        if not rxn:
            raise ValueError("Invalid reaction SMILES")

        # Check reactants, products, and reagents for atom mapping
        for molecule_set in (rxn.GetReactants(), rxn.GetProducts(), rxn.GetAgents()):
            for molecule in molecule_set:
                for atom in molecule.GetAtoms():
                    # Atom mapping is stored in atom.GetAtomMapNum()
                    if atom.GetAtomMapNum() > 0:
                        return True
        return False
    except Exception as e:
        print(f"Error processing reaction SMILES: {e}")
        return False


def reassign_roles(substances: List[dict[str, Any]]) -> List[dict[str, Any]]:
    # First, get the products atom indices
    product_atom_indices: set[int] = set()
    for substance in substances:
        if substance["role"] == "product":

            mol = Chem.MolFromSmiles(substance["smiles"])
            if mol:
                p_indices = set()
                p_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
                product_atom_indices = product_atom_indices.union(p_indices)

    # Now scan reactants and see if any of them should be reagent instead
    for substance in substances:
        if substance["role"] == "reactant":

            mol = Chem.MolFromSmiles(substance["smiles"])
            if mol:
                rc_indices = set()
                rc_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}

                if len(rc_indices.intersection(product_atom_indices)) == 0:
                    substance["role"] = "reagent"

    # Now scan reagents and see if any of them should be reactant instead
    for substance in substances:
        if substance["role"] == "reagent":

            mol = Chem.MolFromSmiles(substance["smiles"])
            if mol:
                rg_indices = set()
                rg_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}

                if len(rg_indices.intersection(product_atom_indices)) > 0:
                    substance["role"] = "reactant"

    return substances


# by ChatGPT
def get_charge(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"[ERROR] SMILES not parseable in the get_charge function. Please check the SMILES {smiles}. Skipping RXSMILES.")
    return Chem.GetFormalCharge(mol)


def reconstruct_rxsmiles(substances: List[dict[str, Any]]) -> str:
    reactants = []
    reagents = []
    products = []

    # Calculate heavy atom count for each substance and assign to the needed array
    for substance in substances:
        mol = Chem.MolFromSmiles(substance["smiles"])
        heavy_atom_count = CalcNumHeavyAtoms(mol)
        substance["heavy_atom_count"] = heavy_atom_count
        if substance["role"] == "reactant":
            reactants.append(substance)
        elif substance["role"] == "reagent":
            reagents.append(substance)
        elif substance["role"] == "product":
            products.append(substance)

    # Sort the substances by heavy atom count so that heavier molecules come first
    reactants_sorted = sorted(reactants, key=lambda x: x["heavy_atom_count"], reverse=True)
    reagents_sorted = sorted(reagents, key=lambda x: x["heavy_atom_count"], reverse=True)
    products_sorted = sorted(products, key=lambda x: x["heavy_atom_count"], reverse=True)

    final_reactants = []
    final_reagents = []
    final_products = []

    final_fragments = []
    total_index = 0
    for substance in reactants_sorted:
        smiles = substance["smiles"]
        final_reactants.append(smiles)
        fragments = smiles.split(".")
        if len(fragments) > 1:
            fragment = ".".join(str(i) for i in range(total_index, total_index + len(fragments)))
            final_fragments.append(fragment)
        total_index += len(fragments)

    for substance in reagents_sorted:
        smiles = substance["smiles"]
        final_reagents.append(smiles)
        fragments = smiles.split(".")
        if len(fragments) > 1:
            fragment = ".".join(str(i) for i in range(total_index, total_index + len(fragments)))
            final_fragments.append(fragment)
        total_index += len(fragments)

    for substance in products_sorted:
        smiles = substance["smiles"]
        final_products.append(smiles)
        fragments = smiles.split(".")
        if len(fragments) > 1:
            fragment = ".".join(str(i) for i in range(total_index, total_index + len(fragments)))
            final_fragments.append(fragment)
        total_index += len(fragments)

    # Reconstruction of the final SMILES
    final_smiles = f"{'.'.join(final_reactants)}>{'.'.join(final_reagents)}>{'.'.join(final_products)}"
    if len(final_fragments) > 0:
        final_smiles += f" |f:{','.join(final_fragments)}|"

    return final_smiles


def normalize_roles(rxsmiles: str) -> str:
    # Throw error if rxsmiles has no atom mapping
    if not rxsmiles_has_atommapping(rxsmiles):
        logger.error(f"[ERROR] No atom mapping found in the RXSMILES: {rxsmiles}.")
        raise RxsmilesAtomMappingException(rxsmiles=rxsmiles, message=f"No atom mapping found in the RXSMILES: {rxsmiles}.")

    # Parse the reaction SMILES into a reaction object
    parsed_rxn = decomposition_utils.parse_reaction_smiles(rxsmiles)

    # Check if the reaction SMILES has a fragment extension
    has_fragment_extension = False
    if parsed_rxn.fragment_groups and len(parsed_rxn.fragment_groups) > 0:
        has_fragment_extension = True
        validate_f_section_indices(rxsmiles)

    rxsmiles_no_extension = parsed_rxn.rxsmiles_no_extension
    fragments_groups = parsed_rxn.fragment_groups
    reactants = parsed_rxn.reactants
    reagents = parsed_rxn.reagents
    products = parsed_rxn.products

    all_substances: List[dict[str, Any]] = (
        [{"smiles": reactant, "role": "reactant", "og_role_index": i} for i, reactant in enumerate(reactants)]
        + [{"smiles": reagent, "role": "reagent", "og_role_index": i} for i, reagent in enumerate(reagents)]
        + [{"smiles": product, "role": "product", "og_role_index": i} for i, product in enumerate(products)]
    )

    for substance in all_substances:
        # Check if substance has more than one fragment
        if "." in substance["smiles"]:
            # order SMILES of fragments so that cations come first and anions last, implementation by ChatGPT
            # This helps to identify the individual reagent substances, despite RDKit does not seem to separate them above the arrow by a plus sign
            # like in the case of reactants and products.
            fragments = substance["smiles"].split(".")
            sorted_fragments = sorted(fragments, key=lambda s: (-get_charge(s), s))
            substance["smiles"] = ".".join(sorted_fragments)

    # Reassign roles based on atom mapping
    reassigned_substances = reassign_roles(all_substances)

    # Reconstruct the reaction SMILES
    reconstructed_rxsmiles = reconstruct_rxsmiles(reassigned_substances)

    return reconstructed_rxsmiles


def compute_rbi(rxsmiles: str) -> float:
    # Parse the reaction SMILES into a reaction object
    parsed_rxn = decomposition_utils.parse_reaction_smiles(rxsmiles)
    return _compute_rbi(parsed_rxn)


def _compute_rbi(parsed_reaction: decomposition_utils.ParsedReaction) -> float:
    rxsmiles = parsed_reaction.rxsmiles

    # Throw error if rxsmiles has no atom mapping
    if not rxsmiles_has_atommapping(parsed_reaction.rxsmiles):
        logger.error(f"[ERROR] No atom mapping found in the RXSMILES: {rxsmiles}.")
        raise RxsmilesAtomMappingException(rxsmiles=rxsmiles, message=f"No atom mapping found in the RXSMILES: {rxsmiles}.")

    all_reactant_atoms: list[Chem.Atom] = []
    all_reactant_atom_indices: set[int] = set()
    reactant_atom_indices: set[int] = set()

    # This is by ChatGPT
    for substance in parsed_reaction.reactants:
        mol = Chem.MolFromSmiles(substance)
        if mol:
            atoms = set([atom for atom in mol.GetAtoms()])
            all_reactant_atoms.extend(atoms)
            reactant_atom_indices = set()
            reactant_atom_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}

            all_reactant_atom_indices = all_reactant_atom_indices.union(reactant_atom_indices)

    # See what fraction of reactant atoms are accounted for in product side
    all_product_atom_indices: set[int] = set()
    product_atom_indices: set[int] = set()

    # This is by ChatGPT
    for substance in parsed_reaction.products:
        mol = Chem.MolFromSmiles(substance)
        if mol:
            product_atom_indices = set()
            product_atom_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
            all_product_atom_indices = all_product_atom_indices.union(product_atom_indices)

    matching_indices = all_product_atom_indices.intersection(all_reactant_atom_indices)

    # Evade div by zero
    if len(all_reactant_atoms) == 0:
        return 0.0

    rbi = float(100.0 * float(len(matching_indices)) / float(len(all_reactant_atoms)))
    return rbi


def compute_pbi(rxsmiles: str) -> float:
    # Parse the reaction SMILES into a reaction object
    parsed_rxn = decomposition_utils.parse_reaction_smiles(rxsmiles)
    return _compute_pbi(parsed_rxn)


def _compute_pbi(parsed_reaction: decomposition_utils.ParsedReaction) -> float:
    rxsmiles = parsed_reaction.rxsmiles

    # Throw error if rxsmiles has no atom mapping
    if not rxsmiles_has_atommapping(parsed_reaction.rxsmiles):
        logger.error(f"[ERROR] No atom mapping found in the RXSMILES: {rxsmiles}.")
        raise RxsmilesAtomMappingException(rxsmiles=rxsmiles, message=f"No atom mapping found in the RXSMILES: {rxsmiles}.")

    all_product_atoms: list[Chem.Atom] = []
    all_product_atom_indices: set[int] = set()
    product_atom_indices: set[int] = set()

    # This is by ChatGPT
    for substance in parsed_reaction.products:
        mol = Chem.MolFromSmiles(substance)
        if mol:
            atoms = set([atom for atom in mol.GetAtoms()])
            all_product_atoms.extend(atoms)
            product_atom_indices = set()
            product_atom_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
            all_product_atom_indices = all_product_atom_indices.union(product_atom_indices)

    # See what fraction of product atoms are accounted for in reactant side
    all_reactant_atom_indices: set[int] = set()
    reactant_atom_indices: set[int] = set()

    # This is by ChatGPT
    for substance in parsed_reaction.reactants:
        mol = Chem.MolFromSmiles(substance)
        if mol:
            reactant_atom_indices = set()
            reactant_atom_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
            all_reactant_atom_indices = all_reactant_atom_indices.union(reactant_atom_indices)

    matching_indices = all_reactant_atom_indices.intersection(all_product_atom_indices)

    # Evade div by zero
    if len(all_product_atoms) == 0:
        return 0.0

    pbi = float(100.0 * float(len(matching_indices)) / float(len(all_product_atoms)))

    return pbi


def compute_tbi(rxsmiles: str) -> float:
    # Parse the reaction SMILES into a reaction object
    parsed_rxn = decomposition_utils.parse_reaction_smiles(rxsmiles)
    return _compute_tbi(parsed_rxn)


def _compute_tbi(parsed_reaction: decomposition_utils.ParsedReaction) -> float:
    rxsmiles = parsed_reaction.rxsmiles

    # Throw error if rxsmiles has no atom mapping
    if not rxsmiles_has_atommapping(parsed_reaction.rxsmiles):
        logger.error(f"[ERROR] No atom mapping found in the RXSMILES: {rxsmiles}.")
        raise RxsmilesAtomMappingException(rxsmiles=rxsmiles, message=f"No atom mapping found in the RXSMILES: {rxsmiles}.")

    # We only need to consider reactant and product atoms, reagents by definition are not accounted for
    all_atoms: list[Chem.Atom] = []
    all_reactant_atom_indices: set[int] = set()
    reactant_atom_indices: set[int] = set()

    # This is by ChatGPT
    for substance in parsed_reaction.reactants:
        mol = Chem.MolFromSmiles(substance)
        if mol:
            atoms = set([atom for atom in mol.GetAtoms()])
            all_atoms.extend(atoms)
            reactant_atom_indices = set()
            reactant_atom_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
            all_reactant_atom_indices = all_reactant_atom_indices.union(reactant_atom_indices)

    # See what fraction of reactant atoms are accounted for in product side
    all_product_atom_indeces: set[int] = set()
    product_atom_indices: set[int] = set()

    # This is by ChatGPT
    for substance in parsed_reaction.products:
        mol = Chem.MolFromSmiles(substance)
        if mol:
            atoms = set([atom for atom in mol.GetAtoms()])
            all_atoms.extend(atoms)
            product_atom_indices = set()
            product_atom_indices = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
            all_product_atom_indeces = all_product_atom_indeces.union(product_atom_indices)

    matching_indices = all_product_atom_indeces.intersection(all_reactant_atom_indices)

    # Evade div by zero
    if len(all_atoms) == 0:
        return 0.0

    # factor of 2 incorporated that matching atom indices are only accounted for only one time up to this point due to the intersection operation
    tbi = float(100.0 * 2.0 * float(len(matching_indices)) / float(len(all_atoms)))
    return tbi


class RxnBalanceOutput(BaseModel):
    rxsmiles: str
    rbi: float
    pbi: float
    tbi: float


def compute_rxn_balance_indices(rxsmiles: str) -> RxnBalanceOutput:
    # Parse the reaction SMILES into a reaction object
    parsed_rxn = decomposition_utils.parse_reaction_smiles(rxsmiles)

    # Compute RBI, PBI, and TBI
    rbi = _compute_rbi(parsed_rxn)
    pbi = _compute_pbi(parsed_rxn)
    tbi = _compute_tbi(parsed_rxn)

    return RxnBalanceOutput(rxsmiles=rxsmiles, rbi=rbi, pbi=pbi, tbi=tbi)
