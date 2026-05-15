#
# The below reaction/molecule drawing code is from ASKCOS. It may of been slightly adapted for AICP purposes.
# Original code can be found here: https://github.com/ncats/aspire_askcos_askcos2_core/blob/main/utils/draw_impl.py
#
import io
import itertools
import logging
import traceback as tb

import rdkit.Chem as Chem
import rdkit.Chem.Draw as Draw
import svgutils.transform
from syngps.utils.decomposition_utils import parse_reaction_smiles
from PIL import Image
from rdkit.Chem import (
    rdAbbreviations,
    rdDepictor,
    rdFMCS,
)
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import rdGeometry

logger = logging.getLogger(__name__)

ABBREVIATIONS = rdAbbreviations.GetDefaultAbbreviations()

HIGHLIGHT_COLORS = [
    (153 / 255, 221 / 255, 255 / 255),  # baby-blue
    (255 / 255, 200 / 255, 153 / 255),  # peach-crayola
    (143 / 255, 237 / 255, 143 / 255),  # light-green
    (255 / 255, 170 / 255, 187 / 255),  # cherry-blossom-pink
    (119 / 255, 170 / 255, 221 / 255),  # iceberg
    (255 / 255, 255 / 255, 170 / 255),  # lemon-yellow-crayola
    (222 / 255, 161 / 255, 222 / 255),  # plum-web
    (153 / 255, 200 / 255, 100 / 255),  # pistachio
]


def get_drawer(w, h, svg=True, transparent=True, options=None):
    """
    Create the appropriate SVG or PNG drawer instance and set options.

    Args:
        w (int): width of drawer in pixels
        h (int): height of drawer in pixels
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        options (MolDrawOptions, optional): custom options object

    Returns:
        MolDraw2DSVG or MolDraw2DCairo instance
    """
    if svg:
        d = rdMolDraw2D.MolDraw2DSVG(w, h)
    else:
        d = rdMolDraw2D.MolDraw2DCairo(w, h)

    if options:
        d.SetDrawOptions(options)

    if transparent:
        d.drawOptions().clearBackground = False

    return d


def get_options(**kwargs):
    """
    Get MolDrawOptions instance with custom defaults.

    Args:
        **kwargs: can be used to set any options value

    Returns:
        rdMolDraw2D.MolDrawOptions
    """
    options = rdMolDraw2D.MolDrawOptions()
    options.additionalAtomLabelPadding = 0.15
    options.fixedBondLength = 16
    options.bondLineWidth = 2
    options.minFontSize = 14
    options.maxFontSize = 14
    for key, val in kwargs.items():
        setattr(options, key, val)
    return options


def get_size(mol):
    """
    Get the size of the input RDKit Mol. Computes 2D coordinates.

    Args:
        mol (Chem.Mol): molecule to get bounds of

    Returns:
        (float, float): x range, y range in molecule units
    """
    if not mol.GetNumConformers():
        rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    xrange = max(xs) - min(xs)
    yrange = max(ys) - min(ys)
    return xrange, yrange


def get_scaled_drawer(mol, svg=True, transparent=True, padding=1.0, options=None):
    """
    Create a drawing canvas which is sized appropriately for the molecule.
    Uses the size of the molecule to get an approximate size for the image.
    RDKit will still reduce the scale if necessary to fit the molecule.

    Args:
        mol (Chem.Mol): molecule object to draw
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        padding (float, optional): extra padding width around border in molecule
            units (default: 1.0)
        options (MolDrawOptions, optional): custom options object

    Returns:
        MolDraw2DSVG or MolDraw2DCairo instance
    """
    xrange, yrange = get_size(mol)

    # Lower bound on size for small molecules
    xrange = max(xrange, 1.5 + 0.75 * xrange) + padding + 2.0
    yrange = max(yrange, 1.5) + padding
    # print(xrange, yrange)

    options = options or get_options()
    scale = options.fixedBondLength
    rd_padding = options.padding
    # print(scale, rd_padding)

    w = int(scale * xrange * (1 + 2 * rd_padding))
    h = int(scale * yrange * (1 + 2 * rd_padding))
    # print(w, h)

    d = get_drawer(w, h, svg=svg, transparent=transparent, options=options)

    return d


def align_molecule(mol, ref):
    """
    Generate conformer that is aligned to reference based on MCS.

    Args:
        mol (Chem.Mol): molecule object to be aligned
        ref (Chem.Mol): reference molecule to align to

    Returns:
        Chem.Mol: Updated mol with aligned conformer
    """
    if not ref.GetNumConformers():
        rdDepictor.Compute2DCoords(ref)

    # Find maximum common substructure between mol and ref
    mcs = rdFMCS.FindMCS(
        [ref, mol],
        matchValences=True,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchChiralTag=True,
    )
    scaffold = Chem.MolFromSmarts(mcs.smartsString)

    # Generate atom mapping for MCS between mol and ref
    mol_match = mol.GetSubstructMatch(scaffold)
    ref_match = ref.GetSubstructMatch(scaffold)
    conf = ref.GetConformer()
    coord_map = {}
    for i, j in zip(ref_match, mol_match):
        pt = conf.GetAtomPosition(i)
        coord_map[j] = rdGeometry.Point2D(pt.x, pt.y)

    # Generate conformer with predetermined coordinates
    rdDepictor.Compute2DCoords(mol, coordMap=coord_map, canonOrient=False, clearConfs=True)
    return mol


def mol_to_image(
    mol,
    svg=True,
    transparent=True,
    padding=1.0,
    return_png=False,
    clear_map=True,
    abbreviate=True,
    bw_atoms=True,
    update=False,
    highlight_atoms=None,
    highlight_bonds=None,
    highlight_atom_colors=None,
    highlight_bond_colors=None,
    show_atom_indices=False,
    reference=None,
    options=None,
    kekulize_mol=False,
    **kwargs,
):
    """
    Create image for the provided RDKit Mol object.

    Args:
        mol (Chem.Mol): molecule object to draw
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        padding (float, optional): extra padding width around border in molecule
            units (default: 1.0)
        return_png (bool, optional): if True, return PNG string instead of PIL.Image
            (default: False)
        clear_map (bool, optional): if True, clear atom map before drawing
            (default: True)
        abbreviate (bool, optional): if True, use functional group abbreviations
            (default: False)
        bw_atoms (bool, optional): if True, use black and white atom palette
            (default: True)
        update (bool, optional): if True, run UpdatePropertyCache before drawing
            (default: False)
        highlight_atoms (list, optional): indices of atoms to highlight
        highlight_bonds (list, optional): indices of bonds to highlight
        highlight_atom_colors (dict, optional): map from atom indices to highlight
            colors
        highlight_bond_colors (dict, optional): map from bond indices to highlight
            colors
        show_atom_indices (bool, optional): if true show atom indices in reaction depiction
        reference (str, optional): reference molecule to align drawing to based on MCS
        options (MolDrawOptions, optional): custom options object
        **kwargs: passed to Draw.PrepareMolForDrawing

    Returns:
        str if svg=True
        PIL.Image instance if svg=False and return_png=False
        bytes if svg=False and return_png=True
    """
    if not mol:
        raise ValueError("Need a valid RDKit molecule to draw!")

    options = options or get_options()

    if bw_atoms:
        options.useBWAtomPalette()
    if show_atom_indices:
        for atom in mol.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetAtomMapNum() if atom.GetAtomMapNum() > 0 else ""))
    if clear_map:
        [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
    if show_atom_indices:
        mol.UpdatePropertyCache(False)
    if abbreviate and not highlight_atoms and not highlight_bonds and clear_map and not show_atom_indices:
        mol.UpdatePropertyCache(False)
        mol = rdAbbreviations.CondenseMolAbbreviations(mol, ABBREVIATIONS)
    if update:
        mol.UpdatePropertyCache(False)  # From legacy code, not sure if truly necessary
    if highlight_atoms or not clear_map:
        options.fixedBondLength = max(options.fixedBondLength, 20)
        padding += 1.0
    if reference:
        reference = Chem.MolFromSmiles(reference)
        if abbreviate and not highlight_atoms and not highlight_bonds and clear_map:
            reference = rdAbbreviations.CondenseMolAbbreviations(reference, ABBREVIATIONS)
        mol = align_molecule(mol, reference)

    # mol = Draw.PrepareMolForDrawing(mol, kekulize=kekulize_mol, **kwargs)
    # Explicit Kekulization handling
    if kekulize_mol:
        work = Chem.Mol(mol)  # copy to avoid side-effects upstream
        try:
            Chem.Kekulize(work, clearAromaticFlags=True)
        except Exception as e:
            logger.warning(f"Kekulize failed ({e}); proceeding without Kekulé form.")
        # Force-clear any residual aromatic flags (belt + suspenders)
        for b in work.GetBonds():
            if b.GetIsAromatic():
                b.SetIsAromatic(False)
        for a in work.GetAtoms():
            if a.GetIsAromatic():
                a.SetIsAromatic(False)
        mol = work
        # Do NOT let PrepareMolForDrawing re-kekulize (it would keep aromatic perception internally)
        mol = Draw.PrepareMolForDrawing(mol, kekulize=False, **kwargs)
    else:
        # Normal path: allow RDKit to keep aromatic perception (classic depiction)
        mol = Draw.PrepareMolForDrawing(mol, kekulize=False, **kwargs)
    arom_bonds = sum(1 for b in mol.GetBonds() if b.GetIsAromatic())
    logger.debug(f"Aromatic bonds after preparation: {arom_bonds} (kekulize_mol={kekulize_mol})")
    if kekulize_mol and arom_bonds > 0:
        logger.debug("Residual aromatic bonds detected; Kekulé depiction may still show aromatic styling.")
    options.prepareMolsBeforeDrawing = False  # Already prepared

    d = get_scaled_drawer(mol, svg=svg, transparent=transparent, padding=padding, options=options)
    if kekulize_mol:
        # Ensure no aromatic circle overlay even if something slipped through
        draw_opts = d.drawOptions()
        if hasattr(draw_opts, "noAromaticCircle"):
            draw_opts.noAromaticCircle = True
    d.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=highlight_bonds,
        highlightAtomColors=highlight_atom_colors,
        highlightBondColors=highlight_bond_colors,
    )
    d.FinishDrawing()

    if svg or return_png:
        return d.GetDrawingText()
    else:
        return Draw._drawerToImage(d)


def draw_arrow(
    arrow_length=80,
    arrow_size=10,
    offset=5,
    padding=10,
    retro=False,
    svg=True,
    transparent=True,
    return_png=False,
):
    """
    Create image of a reaction arrow.

    Args:
        arrow_length (int, optional): length of arrow in pixels (default: 80)
        arrow_size (int, optional): size of arrowhead in pixels (default: 10)
        offset (int, optional): half distance between arrow body lines for retro arrow
            (default: 6)
        padding (int, optional): width of padding around border in pixels (default: 10)
        retro (int, optional): draw arrow style for retro reaction with two body lines
            (default: False)
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        return_png (bool, optional): if True, return PNG string instead of PIL.Image
            (default: False)

    Returns:
        str if svg=True
        PIL.Image instance if svg=False and return_png=False
        bytes if svg=False and return_png=True
    """
    w = arrow_length + 2 * padding
    h = 2 * arrow_size + 2 * padding

    d = get_drawer(w, h, svg=svg, transparent=transparent)
    if not transparent:
        d.ClearDrawing()
        d.SetColour((0, 0, 0))  # Reset drawing color to black

    midh = h // 2
    if retro:
        begin_top = rdGeometry.Point2D(padding, midh + offset)
        end_top = rdGeometry.Point2D(w - padding - offset, midh + offset)
        begin_bottom = rdGeometry.Point2D(padding, midh - offset)
        end_bottom = rdGeometry.Point2D(w - padding - offset, midh - offset)
        top = rdGeometry.Point2D(w - padding - arrow_size, midh + arrow_size)
        bottom = rdGeometry.Point2D(w - padding - arrow_size, midh - arrow_size)
        end = rdGeometry.Point2D(w - padding, midh)
        d.DrawLine(begin_top, end_top, rawCoords=True)
        d.DrawLine(begin_bottom, end_bottom, rawCoords=True)
        d.DrawLine(end, top, rawCoords=True)
        d.DrawLine(end, bottom, rawCoords=True)
    else:
        begin = rdGeometry.Point2D(padding, midh)
        end = rdGeometry.Point2D(w - padding, midh)
        top = rdGeometry.Point2D(w - padding - arrow_size, midh + arrow_size)
        bottom = rdGeometry.Point2D(w - padding - arrow_size, midh - arrow_size)
        d.DrawLine(begin, end, rawCoords=True)
        d.DrawLine(end, top, rawCoords=True)
        d.DrawLine(end, bottom, rawCoords=True)

    d.FinishDrawing()
    if svg or return_png:
        return d.GetDrawingText()
    else:
        return Draw._drawerToImage(d)


def draw_plus(size=20, padding=10, svg=True, transparent=True, return_png=False):
    """
    Create image of a reaction arrow.

    Args:
        size (int, optional): size of plus in pixels (default: 20)
        padding (int, optional): width of padding around border in pixels (default: 10)
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        return_png (bool, optional): if True, return PNG string instead of PIL.Image (default: False)

    Returns:
        str if svg=True
        PIL.Image instance if svg=False and return_png=False
        bytes if svg=False and return_png=True
    """
    w = h = size + 2 * padding

    d = get_drawer(w, h, svg=svg, transparent=transparent)
    if not transparent:
        d.ClearDrawing()
        d.SetColour((0, 0, 0))  # Reset drawing color to black

    left = rdGeometry.Point2D(padding, h // 2)
    right = rdGeometry.Point2D(w - padding, h // 2)
    top = rdGeometry.Point2D(w // 2, h - padding)
    bottom = rdGeometry.Point2D(w // 2, padding)

    d.DrawLine(left, right, rawCoords=True)
    d.DrawLine(top, bottom, rawCoords=True)

    d.FinishDrawing()
    if svg or return_png:
        return d.GetDrawingText()
    else:
        return Draw._drawerToImage(d)


def combine_images_horizontally(images, transparent=True, return_png=False):
    """
    Combines a list of images into a single horizontal image.
    Each image will be vertically centered.
    Images should be provided as SVG strings or PIL.Image objects.

    Args:
        images (list): list of SVG strings or PIL.Image objects to be combined
        transparent (bool, optional): use transparent background (default: True)
        return_png (bool, optional): if True, return PNG string instead of PIL.Image (default: False)

    Returns:
        str if input images are SVG strings
        PIL.Image instance if input images are PIL.Image and return_png=False
        bytes if input images are PIL.Image and return_png=True
    """
    if isinstance(images[0], str):
        svgs = [svgutils.transform.fromstring(s) for s in images]
        heights = [int(svg.height.strip("px")) for svg in svgs]
        height = max(heights)
        widths = [int(svg.width.strip("px")) for svg in svgs]
        width = sum(widths)

        elements = []
        for i, svg in enumerate(svgs):
            offset_x = sum(widths[:i])  # left to right
            offset_y = (height - heights[i]) // 2
            element = svg.getroot()
            element.moveto(offset_x, offset_y)
            elements.append(element)

        # Create a blank image with the full dimensions
        d = get_drawer(width, height, transparent=transparent)
        if not transparent:
            d.ClearDrawing()
        d.FinishDrawing()
        image = svgutils.transform.fromstring(d.GetDrawingText())
        image.append(elements)

        return image.to_str().decode("utf-8")  # Return str instead of bytes
    else:
        heights = [img.size[1] for img in images]
        height = max(heights)
        widths = [img.size[0] for img in images]
        width = sum(widths)
        alpha = 0 if transparent else 255

        # Create a blank image with the full dimensions
        image = Image.new("RGBA", (width, height), (255, 255, 255, alpha))

        for i, img in enumerate(images):
            offset_x = sum(widths[:i])  # left to right
            offset_y = (height - heights[i]) // 2
            image.paste(img, (offset_x, offset_y))

        if return_png:
            bio = io.BytesIO()
            image.save(bio, format="PNG")
            return bio.getvalue()
        else:
            return image


def molecule_smiles_to_image(
    smiles,
    svg=True,
    transparent=True,
    return_png=True,
    highlight=False,
    reacting_atoms=None,
    bonds=False,
    abbreviate=True,
    split=True,
    reference=None,
    show_atom_indices=False,
    kekulize_mol=False,
    **kwargs,
):
    """
    Create image of the provided molecule SMILES string.

    Args:
        smiles (str): SMILES string of molecule to draw
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        return_png (bool, optional): if True, return PNG string instead of PIL.Image (default: False)
        highlight (bool, optional): if True, highlight reacting atoms (default: False)
        reacting_atoms (list, optional): reacting atom info for highlighting
        bonds (bool, optional): highlight bonds between reacting atoms (default: False)
        abbreviate (bool, optional): if True, use functional group abbreviations (default: True)
        split (bool, optional): split and draw molecule fragments separately,
            ignored if highlight=True (default: True)
        reference (str, optional): reference molecule to align drawing to based on MCS
        show_atom_indices (bool, optional): if true show atom indices in reaction depiction
        kekulize_mol (bool, optional): if True, attempt to Kekulize the molecule before drawing (default: True)
        **kwargs: Passed to mol_to_image

    Returns:
        str if svg=True
        PIL.Image instance if svg=False and return_png=False
        bytes if svg=False and return_png=True
    """
    if not highlight and split and smiles.count(".") > 1:
        mols = [Chem.MolFromSmarts(smi) if show_atom_indices else Chem.MolFromSmiles(smi) for smi in smiles.split(".")]

        if kekulize_mol:
            # Kekulize each fragment to avoid drawing issues with aromatic molecules
            for i, mol in enumerate(mols):
                try:
                    Chem.Kekulize(mol, clearAromaticFlags=True)
                    logger.debug("Kekulized molecule fragment successfully.")
                    logger.debug(f"smiles: {Chem.MolToSmiles(mol)}")
                except Exception:
                    raise ValueError(f"Failed to Kekulize molecule: {smiles.split('.')[i]}")

        images = [
            mol_to_image(
                mol,
                svg=svg,
                transparent=transparent,
                abbreviate=abbreviate,
                reference=reference,
                show_atom_indices=show_atom_indices,
                kekulize_mol=kekulize_mol,
                **kwargs,
            )
            for mol in mols
        ]
        return combine_images_horizontally(images, transparent=transparent, return_png=return_png)

    mol = Chem.MolFromSmarts(smiles) if show_atom_indices else Chem.MolFromSmiles(smiles)

    # Kekulize the molecule to avoid drawing issues with aromatic molecules
    if kekulize_mol:
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            raise ValueError(f"Failed to Kekulize molecule: {smiles}")

    # iterate over mol to look for atom index property
    if highlight and reacting_atoms is not None:
        try:
            if len(reacting_atoms) == mol.GetNumAtoms():
                # Reacting atoms data corresponds to site selectivity predictions
                highlight_atoms = list(range(mol.GetNumAtoms()))
                highlight_atom_colors = {x: (1, 1, 1) for x in highlight_atoms}
                options = kwargs.get("options") or get_options()
                scores = reacting_atoms
                for i, atom in enumerate(mol.GetAtoms()):
                    if round(scores[i] * 100) > 5 and atom.GetIsAromatic() and atom.GetTotalNumHs() == 1:
                        highlight_atom_colors[i] = (1 - scores[i], 1, 1 - scores[i])
                        options.atomLabels[i] = "{}%".format(int(round(scores[i] * 100)))
                kwargs["highlight_atoms"] = highlight_atoms
                kwargs["highlight_atom_colors"] = highlight_atom_colors
                kwargs["options"] = options
            else:
                # Reacting atoms data corresponds to atoms to highlight
                atom_idx_map = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms()}
                highlight_atoms = [atom_idx_map[x] for x in reacting_atoms]
                highlight_atom_colors = {x: (0, 1, 0) for x in highlight_atoms}
                kwargs["highlight_atoms"] = highlight_atoms
                kwargs["highlight_atom_colors"] = highlight_atom_colors

                if bonds:
                    # Legacy functionality which previously did not work at all
                    # Updated to work based on guess of original intentions
                    # No effect if reacting atoms are not neighbors
                    highlight_bonds = []
                    for a, b in itertools.combinations(reacting_atoms, 2):
                        bond = mol.GetBondBetweenAtoms(atom_idx_map[a], atom_idx_map[b])
                        if bond is not None:
                            highlight_bonds.append(bond.GetIdx())
                    highlight_bond_colors = {x: (0, 1, 0) for x in highlight_bonds}
                    kwargs["highlight_bonds"] = highlight_bonds
                    kwargs["highlight_bond_colors"] = highlight_bond_colors
        except Exception:
            # Failed to determine highlight colors, continue drawing the molecule without highlighting
            print("Unable to determine atom highlights using specified reacting atoms. Drawing without highlights.")
            tb.print_exc()

    return mol_to_image(
        mol,
        svg=svg,
        transparent=transparent,
        return_png=return_png,
        abbreviate=abbreviate,
        reference=reference,
        show_atom_indices=show_atom_indices,
        kekulize_mol=kekulize_mol,
        **kwargs,
    )


def determine_highlight_colors(mol, frag_map, frag_idx=None):
    """
    Determine highlight colors for reactants and products based on atom map.

    More robust behavior:
      - When frag_idx is provided (building phase): assign every mapped atom's mapno
        to frag_idx in frag_map (overwriting if already present).
      - When frag_idx is None (drawing phase): only highlight atoms whose mapno exists
        in frag_map (skip unknown map numbers instead of raising KeyError).

    Args:
        mol (Chem.Mol): molecule object to analyze
        frag_map (dict[int, int]): mapping from atom map number -> fragment index
        frag_idx (int, optional): current fragment index when building frag_map

    Returns:
        dict: highlight kwargs for mol_to_image
    """
    highlight_atoms = []
    highlight_bonds = []
    highlight_atom_colors = {}
    highlight_bond_colors = {}

    for atom in mol.GetAtoms():
        mapno = atom.GetAtomMapNum()
        if not mapno:
            continue  # RDKit uses 0 for "no mapping"

        # Build-phase: record map number -> fragment index
        if frag_idx is not None:
            frag_map[mapno] = frag_idx

        # Draw-phase: if we don't know this map number, skip it (don't crash)
        frag = frag_map.get(mapno)
        if frag is None:
            continue

        idx = atom.GetIdx()
        color = HIGHLIGHT_COLORS[frag % len(HIGHLIGHT_COLORS)]

        highlight_atoms.append(idx)
        highlight_atom_colors[idx] = color

        # Highlight bonds between same-colored neighboring atoms
        for atom2 in atom.GetNeighbors():
            j = atom2.GetIdx()
            if j < idx and highlight_atom_colors.get(j) == color:
                bond = mol.GetBondBetweenAtoms(idx, j)
                if bond is None:
                    continue
                bond_idx = bond.GetIdx()
                highlight_bonds.append(bond_idx)
                highlight_bond_colors[bond_idx] = color

    return {
        "highlight_atoms": highlight_atoms,
        "highlight_bonds": highlight_bonds,
        "highlight_atom_colors": highlight_atom_colors,
        "highlight_bond_colors": highlight_bond_colors,
    }


def reaction_smiles_to_image(
    smiles,
    svg=True,
    transparent=True,
    return_png=True,
    retro=False,
    highlight=False,
    align=False,
    plus=True,
    update=True,
    show_atom_indices=False,
    kekulize_mols=False,
    **kwargs,
):
    """
    Create image of the provided reaction SMILES string. Omits agents.

    Args:
        smiles (str): SMILES string of reaction to draw
        svg (bool, optional): return SVG image (default: True)
        transparent (bool, optional): use transparent background (default: True)
        return_png (bool, optional): if True, return PNG string instead of PIL.Image (default: False)
        retro (bool, optional): whether this is a retrosynthetic reaction (default: False)
        highlight (bool, optional): if True, highlight based on atom mapping (default: False)
        align (bool, optional): if True, reactants will be aligned to first product (default: False)
        plus (bool, optional): whether include a plus sign between reactants (default: True)
        update (bool, optional): if True, run UpdatePropertyCache before drawing (default: False)
        show_atom_indices (bool, optional): if true show atom indices in reaction depiction
        **kwargs: passed to mol_to_image

    Returns:
        str if svg=True
        PIL.Image instance if svg=False and return_png=False
        bytes if svg=False and return_png=True
    """
    parsed_smiles = parse_reaction_smiles(smiles)
    reactants = parsed_smiles.reactants
    products = parsed_smiles.products
    # agents = parsed_smiles.reagents # Commented out as agents are not currently used
    r_mols = [Chem.MolFromSmiles(r) for r in reactants]
    p_mols = [Chem.MolFromSmiles(p) for p in products]

    if highlight:
        atom_frag_map = {}
        for i, mol in enumerate(p_mols if retro else r_mols):
            determine_highlight_colors(mol, atom_frag_map, i)

    images = []
    for i, mol in enumerate(r_mols):
        if highlight:
            kwargs.update(determine_highlight_colors(mol, atom_frag_map))
        if align:
            kwargs["reference"] = products[0]
        if plus and i > 0:
            images.append(draw_plus(svg=svg, transparent=transparent))
        smiles = Chem.MolToSmarts(mol) if show_atom_indices else Chem.MolToSmiles(mol)
        if (not highlight) and smiles.count(".") > 1:
            images.append(molecule_smiles_to_image(smiles, svg=svg, transparent=transparent, **kwargs))
        else:
            images.append(mol_to_image(mol, svg=svg, transparent=transparent, update=update, show_atom_indices=show_atom_indices, **kwargs))
    images.append(draw_arrow(retro=retro, svg=svg, transparent=transparent))

    for i, mol in enumerate(p_mols):

        if highlight:
            kwargs.update(determine_highlight_colors(mol, atom_frag_map))
        if plus and i > 0:
            images.append(draw_plus(svg=svg, transparent=transparent))
        images.append(
            mol_to_image(mol, svg=svg, transparent=transparent, update=update, show_atom_indices=show_atom_indices, kekulize_mol=kekulize_mols, **kwargs)
        )

    return combine_images_horizontally(images, transparent=transparent, return_png=return_png)
