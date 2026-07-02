import base64

from fastapi import APIRouter, HTTPException
from syngps import (  # moleculeinchi_to_svg,
    is_sdf_parseable,
    molecule_smiles_to_image,
    molecule_smiles_to_svg,
    sdf2smiles,
    smiles2bms,
    smiles2inchikey,
    smiles2sdf,
)
from syngps.errors.errors import InchikeyNotFoundError
from app.app_config import APP_CONFIG
from app.logging_config import logger
from app.models import (
    SimpleSubstanceSmiles2SVGRequest,
    SubstanceInchikey2SVGRequest,
    SubstanceInchikey2SVGResponse,
    SubstanceSDF2SmilesRequest,
    SubstanceSDF2SmilesResponse,
    SubstanceSDFIsValidRequest,
    SubstanceSDFIsValidResponse,
    SubstanceSmiles2BMSRequest,
    SubstanceSmiles2BMSResponse,
    SubstanceSmiles2InchikeyRequest,
    SubstanceSmiles2InchikeyResponse,
    SubstanceSmiles2SDFRequest,
    SubstanceSmiles2SDFResponse,
    SubstanceSmiles2SVGRequest,
    SubstanceSmiles2SVGResponse,
)

##########################
# Substance Utils Router
##########################
substance_utils_router = APIRouter(prefix="/substance_utils", tags=["Substance Utils"])


@substance_utils_router.post("/smiles2inchikey", summary="RDKit SMILES to InChI key")
async def smiles_2_inchikey(request: SubstanceSmiles2InchikeyRequest) -> SubstanceSmiles2InchikeyResponse:
    """
    Retrieves the InChI key of a compound.
    """
    try:
        inchikey = smiles2inchikey(request.smiles)
        return SubstanceSmiles2InchikeyResponse(smiles=request.smiles, inchikey=inchikey)
    except Exception as e:
        logger.error(f"Error converting SMILES to InChI key: {str(e)}")
        raise HTTPException(status_code=400, detail="Error converting SMILES to InChI key")


@substance_utils_router.post("/smiles2bms", summary="RDKit SMILES to Bond Manipulation System")
async def smiles_2_bms(request: SubstanceSmiles2BMSRequest) -> SubstanceSmiles2BMSResponse:
    """
    Retrieves the Bond Manipulation System of a compound.
    """
    try:
        bms = smiles2bms(request.smiles)
        return SubstanceSmiles2BMSResponse(smiles=request.smiles, bms=bms)
    except Exception as e:
        logger.error(f"Error converting SMILES to BMS: {str(e)}")
        raise HTTPException(status_code=500, detail="Error converting SMILES to BMS")


@substance_utils_router.post("/smiles2sdf", summary="RDKit SMILES to Structure Data File")
async def smiles_2_sdf(request: SubstanceSmiles2SDFRequest) -> SubstanceSmiles2SDFResponse:
    """
    Retrieves the Structure Data File of a compound.
    """
    try:
        sdf = smiles2sdf(request.smiles, request.toKekulize)
        return SubstanceSmiles2SDFResponse(smiles=request.smiles, sdf=sdf)
    except Exception as e:
        logger.error(f"Error converting SMILES to SDF: {str(e)}")
        raise HTTPException(status_code=500, detail="Error converting SMILES to SDF")


@substance_utils_router.post("/sdf2smiles", summary="Structure Data File to RDKit SMILES")
async def sdf_2_smiles(request: SubstanceSDF2SmilesRequest) -> SubstanceSDF2SmilesResponse:
    """
    Retrieves the SMILES of a compound from a Structure Data File.
    """
    try:
        smiles = sdf2smiles(request.sdf)
        return SubstanceSDF2SmilesResponse(sdf=request.sdf, smiles=smiles)
    except Exception as e:
        logger.error(f"Error converting SDF to SMILES: {str(e)}")
        raise HTTPException(status_code=500, detail="Error converting SDF to SMILES")


@substance_utils_router.post("/is_sdf_valid", summary="Structure Data File validity check")
async def is_sdf_valid(request: SubstanceSDFIsValidRequest) -> SubstanceSDFIsValidResponse:
    """
    Validates if a Structure Data File is valid for RDKit.
    """
    try:
        valid = is_sdf_parseable(request.sdf)
        return SubstanceSDFIsValidResponse(valid=valid)
    except Exception as e:
        logger.error(f"Error validating SDF: {str(e)}")
        raise HTTPException(status_code=500, detail="Error validating SDF")


#########################
# SVG Utils Router
#########################
svg_utils_router = APIRouter(prefix="/svg")


@svg_utils_router.post("/smiles2svg", summary="Generate SVG for a given substance smiles")
async def smiles_2_svg(request: SubstanceSmiles2SVGRequest) -> SubstanceSmiles2SVGResponse:
    """
    Generates SVG for a given substance smiles.

    **Args:**
    - `smiles` (str): The molecule smiles to be converted to SVG.
    - `monochrome_atoms` (bool): If True, the atoms will be black.
    - `width` (int): The width of the SVG image.
    - `height` (int): The height of the SVG image.
    - `base64_encode` (bool): If True, the SVG will be base64 encoded.

    **Returns:**
    - `SubstanceSmiles2SVGResponse`: The SVG representation of the molecule. Contains 'smiles' and 'svg' fields.

    """
    try:
        svg = moleculesmiles_to_svg(request.smiles, monochrome_atoms=request.monochrome_atoms, img_width=request.width, img_height=request.height)
        svg = svg.replace('"', "'")
        if request.base64_encode:
            svg = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
        return SubstanceSmiles2SVGResponse(smiles=request.smiles, svg=svg, base64_encoded=request.base64_encode)
    except Exception as e:
        logger.error(f"Error converting SMILES to SVG: {str(e)}")
        raise HTTPException(status_code=500, detail="Error converting SMILES to SVG") from e


@svg_utils_router.post("/simple_smiles2svg", summary="Generate SVG for a given substance smiles using ASKCOS method")
async def simple_smiles_2_svg(request: SimpleSubstanceSmiles2SVGRequest) -> SubstanceSmiles2SVGResponse:
    """
    Generates SVG for a given substance smiles.

    **Args:**
    - `smiles` (str): The molecule smiles to be converted to SVG.
    - `base64_encode` (bool): If True, the SVG will be base64 encoded.

    **Returns:**
    - `SubstanceSmiles2SVGResponse`: The SVG representation of the molecule. Contains 'smiles' and 'svg' fields.

    """
    try:
        svg = molecule_smiles_to_image(request.smiles, transparent=False, kekulize_mol=request.kekulize_mol, show_atom_indices=False)
        svg = svg.replace('"', "'")
        if request.base64_encode:
            svg = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
        return SubstanceSmiles2SVGResponse(smiles=request.smiles, svg=svg, base64_encoded=request.base64_encode)
    except Exception as e:
        logger.error(f"Error converting SMILES to SVG: {str(e)}")
        raise HTTPException(status_code=500, detail="Error converting SMILES to SVG") from e


# Include the routers in the main router
substance_utils_router.include_router(svg_utils_router)
