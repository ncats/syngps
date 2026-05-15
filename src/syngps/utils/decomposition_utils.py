# Author: Gergely Zahoranszky-Kohalmi PhD, Nathan Miller
#
# Email: gergely.zahoranszky-kohalmi@nih.gov, nathan.miller@nih.gov
#
# Organization: NCATS/NIH
#
# Aim: Parses a reaction SMILES string into its components (reactants, reagents, products) based on fragment groups.
#
#
import re
from typing import List, Optional, Tuple

from pydantic import BaseModel


class RxSmilesInput(BaseModel):
    rxsmiles: str


class Component(BaseModel):
    index: int
    molecule: str


class ReactionComponent(BaseModel):
    smiles: str
    original_index: int


class ReactionComponents(BaseModel):
    reactants: List[ReactionComponent]
    reagents: List[ReactionComponent]
    products: List[ReactionComponent]


class FragmentGroup(BaseModel):
    groups: List[List[int]]


class ParsedReaction(BaseModel):
    rxsmiles: str
    fragment_groups: Optional[List[List[int]]]
    reactants: List[str]
    reagents: List[str]
    products: List[str]

    @property
    def reactant_smiles(self):
        return ".".join(self.reactants)

    @property
    def reagent_smiles(self):
        return ".".join(self.reagents)

    @property
    def product_smiles(self):
        return ".".join(self.products)

    @property
    def rxsmiles_no_extension(self):
        return self.rxsmiles.split("|")[0].strip()


class FragmentGroupError(Exception):
    """Custom exception for invalid fragment groups crossing component boundaries."""

    pass


def has_extension(rxsmiles: str) -> bool:
    """
    Check if the reaction SMILES string has an extension.

    Args:
        rxsmiles (str): The reaction SMILES string.

    Returns:
        bool: True if the SMILES string has an extension, False otherwise.
    """
    return "|" in rxsmiles


def has_fragment_extension(rxsmiles: str) -> bool:
    """
    Check if the reaction SMILES string has a fragment extension.

    Args:
        rxsmiles (str): The reaction SMILES string.

    Returns:
        bool: True if the SMILES string has a fragment extension, False otherwise.
    """
    return "f:" in rxsmiles


def parse_fragment_groups(extension: str) -> FragmentGroup:
    """
    Parse fragment groups from the extension string.

    Args:
        extension (str): The extension string containing fragment group information.

    Returns:
        FragmentGroup: Parsed groups of fragment indices.
    """
    # Updated regex to stop at `|` or `^` and handle trailing commas
    fragment_match = re.search(r"f:([^|^,]+(?:,[^|^,]+)*)", extension)
    if fragment_match:
        # Process the `f:` portion
        fragment_groups = fragment_match.group(1).strip().split(",")
        parsed_groups = [[int(index) for index in group.split(".")] for group in fragment_groups if group]
        return FragmentGroup(groups=parsed_groups)
    else:
        return FragmentGroup(groups=[])


def extract_smiles(rxsmiles: str) -> str:
    """
    Extract the SMILES part from the reaction SMILES string.

    Args:
        rxsmiles (str): The reaction SMILES string.

    Returns:
        str: The SMILES part of the reaction string.
    """
    return rxsmiles.split("|")[0].strip()


def extract_extension(rxsmiles: str) -> str:
    """
    Extract the extension part from the reaction SMILES string.

    Args:
        rxsmiles (str): The reaction SMILES string.

    Returns:
        str: The extension part of the reaction string.
    """
    return "|" + rxsmiles.split("|")[1].strip()


def extract_smiles_and_extension(rxsmiles: str) -> Tuple[str, str]:
    """
    Extract SMILES and the fragment extension.

    Args:
        rxsmiles (str): The reaction SMILES string with an optional extension.

    Returns:
        Tuple[str, str]: A tuple containing the SMILES part and the extension part.
    """
    parts = rxsmiles.split("|", 1)
    smiles = parts[0].strip()
    extension = "|" + parts[1].strip() if len(parts) > 1 else ""
    return smiles, extension


def extract_original_reaction_components(smiles: str) -> ReactionComponents:
    """
    Extract reactants, reagents, and products from the SMILES string. Keep them as separate lists,
    each including the original index of the component.

    Args:
        smiles (str): The reaction SMILES string without any extensions.

    Returns:
        ReactionComponents: Parsed reactants, reagents, and products.
    """
    # Ensure smiles has atleast 2 '>'
    if smiles.count(">") < 2:
        raise ValueError("Invalid reaction SMILES string. It should contain at least 2 '>' characters.")

    # Split the SMILES string into reactants, reagents, and products
    components = smiles.split(">")

    reactants = (
        [ReactionComponent(smiles=smile.strip(), original_index=index) for index, smile in enumerate(components[0].strip().split(".")) if components[0]]
        if len(components) > 0
        else []
    )

    if len(reactants) == 0:
        raise ValueError("Invalid reaction SMILES string. It should contain at least one reactant.")

    reagents = (
        [
            ReactionComponent(smiles=smile.strip(), original_index=index + len(reactants))
            for index, smile in enumerate(components[1].strip().split("."))
            if components[1]
        ]
        if len(components) > 1
        else []
    )

    products = (
        [
            ReactionComponent(smiles=smile.strip(), original_index=index + len(reactants) + len(reagents))
            for index, smile in enumerate(components[2].strip().split("."))
            if components[2]
        ]
        if len(components) > 2
        else []
    )

    if len(products) == 0:
        raise ValueError("Invalid reaction SMILES string. It should contain at least one product.")

    return ReactionComponents(reactants=reactants, reagents=reagents, products=products)


def validate_fragment_groups(fragment_groups: List[List[int]], reaction_components: ReactionComponents) -> None:
    """
    Validate that each fragment group does not cross component boundaries.

    Args:
        fragment_groups (List[List[int]]): List of fragment groups, where each group is a list of indices.
        reaction_components (ReactionComponents): Parsed reactants, reagents, and products.

    Raises:
        FragmentGroupError: If any group spans across different components or if all booleans are false.
    """
    reactant_indices = [reactant.original_index for reactant in reaction_components.reactants]
    reagent_indices = [reagent.original_index for reagent in reaction_components.reagents]
    product_indices = [product.original_index for product in reaction_components.products]

    for group in fragment_groups:
        within_reactants = all(i in reactant_indices for i in group)
        within_reagents = all(i in reagent_indices for i in group)
        within_products = all(i in product_indices for i in group)

        # Check if exactly one of the three is True
        true_count = sum([within_reactants, within_reagents, within_products])

        if true_count != 1:
            raise FragmentGroupError("Invalid fragment group. Fragments either span across components or are invalid.")


def join_reaction_components(reaction_components: ReactionComponents, fragment_groups: List[List[int]]) -> Tuple[List[str], List[str], List[str]]:
    """
    Join reactants, reagents, and products by '.' according to the fragment groups,
    and retain any reaction components not part of a fragment group.

    Args:
        reaction_components (ReactionComponents): Parsed reactants, reagents, and products.
        fragment_groups (List[List[int]]): List of fragment groups, where each group is a list of indices.

    Returns:
        ParsedReaction: An object containing lists of joined reactants, reagents, and products.
    """
    # Create a dictionary to map original indices to their SMILES
    index_to_smiles = {comp.original_index: comp.smiles for comp in reaction_components.reactants + reaction_components.reagents + reaction_components.products}

    # Track indices that are part of a fragment group
    grouped_indices = set(i for group in fragment_groups for i in group)

    # Initialize lists to collect grouped SMILES for reactants, reagents, and products
    grouped_reactants = []
    grouped_reagents = []
    grouped_products = []

    # Process each fragment group
    for group in fragment_groups:
        try:
            # Map each index in the group to its SMILES and join by '.'
            group_smiles = ".".join(index_to_smiles[i] for i in group)
        except KeyError as e:
            raise FragmentGroupError(f"Index {e.args[0]} not found in reaction components.") from e

        # Determine which component list this group belongs to
        if all(i in [r.original_index for r in reaction_components.reactants] for i in group):
            grouped_reactants.append(group_smiles)
        elif all(i in [re.original_index for re in reaction_components.reagents] for i in group):
            grouped_reagents.append(group_smiles)
        elif all(i in [p.original_index for p in reaction_components.products] for i in group):
            grouped_products.append(group_smiles)
        else:
            raise FragmentGroupError("Fragment groups cannot span across reactants, reagents, and products.")

    # Add ungrouped components to their respective lists
    ungrouped_reactants = [r.smiles for r in reaction_components.reactants if r.original_index not in grouped_indices]
    ungrouped_reagents = [re.smiles for re in reaction_components.reagents if re.original_index not in grouped_indices]
    ungrouped_products = [p.smiles for p in reaction_components.products if p.original_index not in grouped_indices]

    # Combine grouped and ungrouped components
    final_reactants = grouped_reactants + ungrouped_reactants
    final_reagents = grouped_reagents + ungrouped_reagents
    final_products = grouped_products + ungrouped_products

    # Return the ParsedReaction object with grouped reactants, reagents, and products
    return (final_reactants, final_reagents, final_products)


def parse_reaction_smiles(rxsmiles: str) -> ParsedReaction:
    """
    Parse the reaction SMILES string into components based on fragments.

    Args:
        rxsmiles (str): The reaction SMILES string with an optional extension.

    Returns:
        ParsedReaction: Parsed reaction components (reactants, reagents, products).
    """
    # Step 1: Extract the SMILES and extension
    smiles, extension = extract_smiles_and_extension(rxsmiles)

    # Step 2: Extract reaction components
    components = extract_original_reaction_components(smiles)

    # Step 2b: If there are no fragment groups, return the components as is
    if not extension:
        return ParsedReaction(
            rxsmiles=rxsmiles,
            fragment_groups=None,
            reactants=[comp.smiles for comp in components.reactants],
            reagents=[comp.smiles for comp in components.reagents],
            products=[comp.smiles for comp in components.products],
        )

    # Step 3: Parse fragment groups
    fragment_groups = parse_fragment_groups(extension).groups

    # Step 4: Validate fragment groups to ensure they do not span across reactants, reagents, and products
    validate_fragment_groups(fragment_groups, components)

    # Step 5: Join reaction components based on fragment groups
    reactants, reagents, products = join_reaction_components(components, fragment_groups)

    # Step 6: Create the final reaction SMILES
    return ParsedReaction(rxsmiles=rxsmiles, fragment_groups=fragment_groups, reactants=reactants, reagents=reagents, products=products)
