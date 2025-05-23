from collections import Counter
from enum import Enum
from typing import List

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Mol
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from errors import NotCanonicalizableSmilesException
import logging
logger = logging.getLogger(__name__)


class RxnSvgDepictionMode(Enum):
    SIMPLE = "simple"
    ATOM_MAP = "atom_map"
    HIGHLIGHT_WO_INDICES = "highlight_wo_indices"
    HIGHLIGHT_WITH_INDICES = "highlight_with_indices"


def smiles2inchikey(smiles: str) -> str:
    """
    Convert SMILES string input to InChIKey string.

    Args:
        smiles (str): SMILES string.

    Returns:
        str: InChIKey string.

    Raises:
        ValueError: If conversion fails or input is invalid.
    """
    if not smiles or smiles == "":
        raise ValueError("Empty SMILES string provided")

    try:
        mol = Chem.MolFromSmiles(smiles)
        inchikey = Chem.inchi.MolToInchiKey(mol)
    except Exception:
        logger.error(f"Failed to convert SMILES to InChIKey: {smiles}", exc_info=True)
        raise ValueError(f"Failed to convert SMILES to InChIKey: {smiles}")
    return inchikey


def smiles2bms(smiles: str) -> str:
    """
    Convert SMILES string input to BMS (Bemis-Murcko scaffold) SMILES string.

    Args:
        smiles (str): SMILES string.

    Returns:
        str: BMS SMILES string.

    Raises:
        ValueError: If conversion fails or input is invalid.
    """
    if not smiles or smiles == "":
        raise ValueError("Empty SMILES string provided")

    try:
        mol = Chem.MolFromSmiles(smiles)
        bms_smiles = Chem.MolToSmiles(GetScaffoldForMol(mol))
    except Exception:
        logger.error(f"Failed to convert SMILES to BMS: {smiles}", exc_info=True)
        raise ValueError(f"Failed to convert SMILES to BMS: {smiles}")
    return bms_smiles


def smiles2sdf(smiles: str, toKekulize: bool = False) -> str:
    """
    Convert SMILES string input to SDF (Structure Data File) string.

    Args:
        smiles (str): SMILES string.
        toKekulize (bool, optional): Kekulize the molecule. Defaults to False.

    Returns:
        str: SDF string

    Raises:
        ValueError: If conversion fails or input is invalid.
    """
    if not smiles or smiles == "":
        raise ValueError("Empty SMILES string provided")

    try:
        mol = Chem.MolFromSmiles(smiles)
        sdf = Chem.MolToMolBlock(mol, kekulize=toKekulize)
    except Exception:
        logger.error(f"Failed to convert SMILES to SDF: {smiles}", exc_info=True)
        raise ValueError(f"Failed to convert SMILES to SDF: {smiles}")
    return sdf


def sdf2smiles(sdf: str) -> str:
    """
    Convert SDF (Structure Data File) string input to SMILES string.

    Args:
        sdf: SDF string

    Returns:
        smiles (str): SMILES string.

    Raises:
        ValueError: If conversion fails or input is invalid.
    """
    if not sdf or sdf == "":
        raise ValueError("Empty SDF provided")

    try:
        mol = Chem.MolFromMolBlock(sdf)
        smiles = Chem.MolToSmiles(mol)
    except Exception:
        logger.error(f"Failed to convert SDF to SMILES: {sdf}", exc_info=True)
        raise ValueError(f"Failed to convert SDF to SMILES: {sdf}")
    return smiles


def is_sdf_parseable(sdf: str) -> bool:
    """
    Checks whether the input SDF is parseable by RDKit.

    Args:
        sdf: SDF string

    Returns:
        bool: True if valid, False otherwise

    Raises:
        ValueError: If the input is erroneous
    """
    if not sdf or sdf == "":
        raise ValueError("Empty SDF provided")

    try:
        mol = Chem.MolFromMolBlock(sdf)
    except Exception:
        logger.error(f"Failed to validate SDF: {sdf}", exc_info=True)
        raise ValueError(f"Failed to validate SDF: {sdf}")
    return mol is not None


def parse_rxn_extended_smiles(extended_rxn_smiles: str) -> AllChem.ChemicalReaction:
    """
    Parse reaction SMILES string.

    Args:
        rxn_smiles (str): Reaction SMILES string.

    Returns:
        rdChemReactions.ChemicalReaction: RDKit ChemicalReaction object.

    Raises:
        ValueError: If conversion fails.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(extended_rxn_smiles)
        rxn.Initialize()
    except Exception:
        logger.error(f"Failed to parse reaction: {extended_rxn_smiles}", exc_info=True)
        raise ValueError(f"Failed to parse reaction: {extended_rxn_smiles}")
    return rxn


def is_reaction_valid(extended_rxn_smiles: str) -> bool:
    """
    Check if reaction SMILES is valid.

    Args:
        rxn_smiles (str): Reaction SMILES string.

    Returns:
        bool: True if valid, False otherwise.

    Raises:
        ValueError: If conversion fails.
    """
    rxn = parse_rxn_extended_smiles(extended_rxn_smiles)
    (warning_count, error_count) = rxn.Validate()

    if warning_count > 0:
        logger.debug(f"{warning_count} warnings found for reaction: {extended_rxn_smiles}")
    if error_count > 0:
        logger.error(f"{error_count} errors found for reaction: {extended_rxn_smiles}")
        return False
    return True


def is_reaction_balanced(extended_rxn_smiles: str) -> bool:
    """
    Check if reaction SMILES is balanced, with improved handling for common molecules.

    Args:
        extended_rxn_smiles (str): Reaction SMILES string.

    Returns:
        bool: True if balanced, False otherwise.
    """
    # Assume this function exists and returns a boolean
    is_rxn_valid = is_reaction_valid(extended_rxn_smiles)

    if not is_rxn_valid:
        # Assume logger is defined elsewhere
        logger.debug(f"Skipping balance check for invalid reaction: {extended_rxn_smiles}")
        return False

    # Assume this function exists and returns an RDKit reaction object
    rxn = parse_rxn_extended_smiles(extended_rxn_smiles)

    reactants_atoms: Counter = Counter()
    products_atoms: Counter = Counter()

    # Helper function to count atoms, ignoring hydrogens
    def count_atoms(mols, atoms_counter):
        for mol in mols:
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 1:  # Ignore hydrogen
                    atoms_counter[atom.GetSymbol()] += 1

    count_atoms(rxn.GetReactants(), reactants_atoms)
    count_atoms(rxn.GetProducts(), products_atoms)

    return reactants_atoms == products_atoms


def rxsmiles_to_svg(
    rxsmiles: str,
    depiction_mode: RxnSvgDepictionMode = RxnSvgDepictionMode.SIMPLE,
    monochrome_atoms: bool = False,
    img_width: int = 1800,
    img_height: int = 600,
) -> str:
    """
    Convert a reaction SMILES string into an SVG representation. Depiction mode can be set to SIMPLE, ATOM_MAP, HIGHLIGHT_WO_INDICES, or HIGHLIGHT_WITH_INDICES.

    Args:
        rxsmiles (str): The SMILES string representation of the reaction.
        depiction_mode (RxnSvgDepictionMode, optional): Depiction mode for the reaction. Defaults to RxnSvgDepictionMode.SIMPLE.
        monochrome_atoms (bool, optional): Flag to render atoms in monochrome. Defaults to False.
        img_width (int, optional): Width of the resulting image. Defaults to 1800.
        img_height (int, optional): Height of the resulting image. Defaults to 600.

    Returns:
        str: SVG string representation of the reaction.

    Raises:
        ValueError: If the reaction cannot be parsed from the SMILES string.
    """
    rxn = AllChem.ReactionFromSmarts(rxsmiles, useSmiles=True)

    # Set atom map numbers if needed
    if depiction_mode == RxnSvgDepictionMode.SIMPLE:
        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)  # Set atom map number starting from 1
        for mol in rxn.GetProducts():
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)  # Set atom map number starting from 1

    elif depiction_mode == RxnSvgDepictionMode.HIGHLIGHT_WITH_INDICES:
        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetAtomMapNum()))  # Set atom map number starting from 1
        for mol in rxn.GetProducts():
            for atom in mol.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetAtomMapNum()))  # Set atom map number starting from 1

    # Initialize SVG drawer
    d2d = Draw.MolDraw2DSVG(img_width, img_height)

    # Set monochrome atoms
    if monochrome_atoms:
        d2d.drawOptions().updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict})

    # Draw reaction
    if depiction_mode in [RxnSvgDepictionMode.HIGHLIGHT_WO_INDICES, RxnSvgDepictionMode.HIGHLIGHT_WITH_INDICES]:
        # Highlight fragment integration
        d2d.DrawReaction(rxn, highlightByReactant=True)
    else:
        # Show atom mapping
        d2d.DrawReaction(rxn)

    # Finish drawing
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    return svg


def moleculesmiles_to_svg(mol_smiles: str, monochrome_atoms: bool = False, kekulize_mol: bool = False, img_width: int = 1800, img_height: int = 600) -> str:
    """
    Convert a molecule SMILES string into an SVG representation.

    Args:
        mol_smiles (str): The SMILES string representation of the molecule.
        img_width (int, optional): Width of the resulting image. Defaults to 1800.
        img_height (int, optional): Height of the resulting image. Defaults to 600.

    Returns:
        str: SVG string representation of the molecule.

    Raises:
        ValueError: If the molecule cannot be parsed from the SMILES string.
    """
    # Check if SMILES string is empty
    if mol_smiles == "":
        logger.error("Empty SMILES string provided")
        raise ValueError("Empty SMILES string provided")

    # Create molecule from SMILES
    mol = Chem.MolFromSmiles(mol_smiles)

    # Kekulize molecule if needed
    if kekulize_mol:
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Chem.rdchem.KekulizeException:
            logger.error(f"Failed to kekulize molecule: {mol_smiles}", exc_info=True)
            raise ValueError(f"Failed to kekulize molecule: {mol_smiles}")

    # Check if molecule is valid
    if mol is None:
        logger.error(f"Invalid SMILES string: {mol_smiles}")
        raise ValueError(f"Invalid SMILES string: {mol_smiles}")

    # Initialize SVG drawer
    d2d = Draw.MolDraw2DSVG(img_width, img_height)

    # Set monochrome atoms
    if monochrome_atoms:
        d2d.drawOptions().updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict})

    # Draw molecule
    try:
        d2d.DrawMolecule(mol)
    except Exception as e:
        logger.error(f"Failed to draw molecule: {mol_smiles}", exc_info=True)
        raise ValueError(f"Failed to draw molecule: {mol_smiles} due to {str(e)}") from e
    d2d.FinishDrawing()

    # Return SVG
    svg = d2d.GetDrawingText()
    return svg


def moleculeinchi_to_svg(mol_inchi: str, monochrome_atoms: bool = False, kekulize_mol: bool = False, img_width: int = 1800, img_height: int = 600) -> str:
    """
    Convert a molecule InChI string into an SVG representation.

    Args:
        mol_inchi (str): The InChI string representation of the molecule.
        img_width (int, optional): Width of the resulting image. Defaults to 1800.
        img_height (int, optional): Height of the resulting image. Defaults to 600.

    Returns:
        str: SVG string representation of the molecule.

    Raises:
        ValueError: If the molecule cannot be parsed from the InChI string.
    """
    # Check if InChI string is empty
    if not mol_inchi:
        logger.error("Empty InChI string provided")
        raise ValueError("Empty InChI string provided")

    # Create molecule from InChI
    mol = Chem.MolFromInchi(mol_inchi)

    # Kekulize molecule if needed
    if kekulize_mol:
        try:
            Chem.Kekulize(mol)
        except Chem.rdchem.KekulizeException:
            logger.error(f"Failed to kekulize molecule: {mol_inchi}", exc_info=True)
            raise ValueError(f"Failed to kekulize molecule: {mol_inchi}")

    # Check if molecule is valid
    if mol is None:
        logger.error(f"Invalid InChI string: {mol_inchi}")
        raise ValueError(f"Invalid InChI string: {mol_inchi}")

    # Initialize SVG drawer
    d2d = Draw.MolDraw2DSVG(img_width, img_height)

    # Set monochrome atoms
    if monochrome_atoms:
        d2d.drawOptions().updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict})

    # Draw molecule
    try:
        d2d.DrawMolecule(mol)
    except Exception as e:
        logger.error(f"Failed to draw molecule: {mol_inchi}", exc_info=True)
        raise ValueError(f"Failed to draw molecule: {mol_inchi} due to {str(e)}") from e
    d2d.FinishDrawing()

    # Return SVG
    svg = d2d.GetDrawingText()
    return svg


def canonicalize_molecule_smiles(molecule_smiles: str, remove_atom_mapping: bool = False) -> str:
    """
    Canonicalize a SMILES string to its standard form.

    Args:
        molecule_smiles (str): The SMILES string to canonicalize.
        remove_atom_mapping (bool): Flag to remove atom-mapping numbers from the SMILES string.

    Returns:
        str: The canonicalized SMILES string.

    Raises:
        NotCanonicalizableSmilesException: If the SMILES string cannot be converted to a molecule.
    """
    mol: Mol = Chem.MolFromSmiles(molecule_smiles)
    if not mol:
        raise NotCanonicalizableSmilesException(molecule_smiles)
    if remove_atom_mapping:
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                atom.ClearProp("molAtomMapNumber")
    return Chem.MolToSmiles(mol, canonical=True)


def process_reaction(rxn: str) -> str:
    """
    Process and canonicalize a chemical reaction represented in SMILES format.

    Args:
        rxn (str): The reaction SMILES string, formatted as 'reactants>reagents>products'.

    Returns:
        str: The canonicalized reaction SMILES string or an empty string if any component cannot be canonicalized.

    Raises:
        ValueError: If the reaction cannot be canonicalized.
    """
    reactants, reagents, products = rxn.split(">")
    try:
        precursors: List[str] = [canonicalize_molecule_smiles(r, True) for r in reactants.split(".")]
        if reagents.strip():
            precursors += [canonicalize_molecule_smiles(r, True) for r in reagents.split(".")]
        product_list: List[str] = [canonicalize_molecule_smiles(p, True) for p in products.split(".")]
    except NotCanonicalizableSmilesException as e:
        logger.error(f"Failed to canonicalize reaction: {rxn}", exc_info=True)
        raise ValueError(f"Failed to canonicalize reaction: {rxn}") from e

    joined_precursors = ".".join(sorted(precursors))
    joined_products = ".".join(sorted(product_list))
    return f"{joined_precursors}>>{joined_products}"
