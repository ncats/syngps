import base64

from fastapi import APIRouter, HTTPException
from syngps import (
    is_reaction_balanced,
    is_reaction_valid,
    reaction_smiles_to_image,
    rxsmiles_to_svg,
)
from syngps.errors import RxsmilesAtomMappingException
from syngps.utils import (
    compute_rxn_balance_indices,
    normalize_roles,
    parse_reaction_smiles,
)
from app.app_config import APP_CONFIG
from app.logging_config import logger
from app.models import (
    AtommapRequest,
    AtommapResponse,
    ComputeBalanceIndicesRequest,
    ComputeBalanceIndicesResponse,
    HzAtommapRequest,
    HzAtommapResponse,
    HzNormalizeRoleRequest,
    HzNormalizeRoleResponse,
    HzRxnoRequest,
    HzRxnoResponse,
    NormalizeRoleRequest,
    NormalizeRoleResponse,
    ParseRxsmilesResponse,
    RxnClassifyRequest,
    RxnClassifyResponse,
    RxnIsBalancedRequest,
    RxnIsBalancedResponse,
    RxnIsValidRequest,
    RxnIsValidResponse,
    Rxsmiles2SVGRequest,
    Rxsmiles2SVGResponse,
    RxsmilesRequest,
    SimpleRxsmiles2SVGRequest,
    YieldPredictionRequest,
    YieldPredictionResponse,
)
from app.routers.asckos_wrappers import get_atom_mapping_rxnmapper


def _hazelnut_not_available_http_exc() -> HTTPException:
    return HTTPException(
        status_code=501,
        detail=(
            "This endpoint requires NextMove Hazelnut (licensed software) and is not available in this open-source repository. "
            "Use a licensed Hazelnut deployment to enable this functionality."
        ),
    )

#########################
# Reaction Utils Router
#########################
reaction_utils_router = APIRouter(prefix="/reaction_utils", tags=["Reaction Utils"])


@reaction_utils_router.post("/parse_rxnsmiles", summary="Parse RXN Smiles into its components")
async def parse_rxnsmiles(request: RxsmilesRequest) -> ParseRxsmilesResponse:
    """
    Parses a RXN Smiles string into its components.
    """
    try:
        parsed_rxn = parse_reaction_smiles(request.rxsmiles)

        return ParseRxsmilesResponse(
            rxsmiles=request.rxsmiles,
            reactants=parsed_rxn.reactants,
            reagents=parsed_rxn.reagents,
            products=parsed_rxn.products,
            fragment_groups=parsed_rxn.fragment_groups,
        )
    except Exception as e:
        logger.error(f"Error parsing RXN Smiles: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal error parsing RXN Smiles")


@reaction_utils_router.post("/normalize_roles", summary="Normalize reaction roles from a RXN Smiles")
async def normalize_rxsmiles_roles(request: NormalizeRoleRequest) -> NormalizeRoleResponse:
    """
    Normalizes the roles of a reaction from a RXN Smiles string. Input string must be a valid RXN Smiles
    with atom mapping.
    """
    rxsmiles = request.rxsmiles
    if request.atom_map:
        atom_map_response = await get_atom_map(AtommapRequest(smiles=rxsmiles))
        rxsmiles = atom_map_response.mapped_rxsmiles

    try:
        normalized_rxn = normalize_roles(rxsmiles)
        return NormalizeRoleResponse(original_rxsmiles=request.rxsmiles, rxsmiles=normalized_rxn, atom_mapped=request.atom_map)
    except RxsmilesAtomMappingException:
        raise HTTPException(status_code=400, detail="Error parsing RXN Smiles: Atom mapping required")
    except Exception as e:
        logger.error(f"Error normalizing roles: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal error normalizing roles")


@reaction_utils_router.post("/normalize_and_parse", summary="Normalize reaction roles from a RXN Smiles and then parses the RXN Smiles")
async def normalize_and_parse_rxsmiles(request: NormalizeRoleRequest) -> ParseRxsmilesResponse:
    """
    Normalizes the roles of a reaction from a RXN Smiles string and then parses the RXN Smiles.
    """
    normalized_rxn = await normalize_rxsmiles_roles(request)
    try:
        parsed_rxn = parse_reaction_smiles(normalized_rxn.rxsmiles)

        return ParseRxsmilesResponse(
            rxsmiles=normalized_rxn.rxsmiles,
            reactants=parsed_rxn.reactants,
            reagents=parsed_rxn.reagents,
            products=parsed_rxn.products,
            fragment_groups=parsed_rxn.fragment_groups,
        )
    except Exception as e:
        logger.error(f"Error normalizing and parsing RXN Smiles: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal error normalizing and parsing RXN Smiles")


@reaction_utils_router.post("/compute_balance_indices", summary="Calculates the balance indices of a reaction. NEED REFERENCES FROM GERGELY.")
async def compute_balance_indices(request: ComputeBalanceIndicesRequest) -> ComputeBalanceIndicesResponse:
    """
    Calculates the balance indices of a reaction.
    """
    rxsmiles = request.rxsmiles
    if request.atom_map:
        try:
            atom_map_response = await get_atom_map(AtommapRequest(smiles=rxsmiles))
            rxsmiles = atom_map_response.mapped_rxsmiles
        except Exception:
            logger.error("Error calculating balance indices.", exc_info=True)
            raise HTTPException(status_code=400, detail="Error computing balance indices")

    try:
        balance_indices = compute_rxn_balance_indices(rxsmiles)
        return ComputeBalanceIndicesResponse(
            rxsmiles=rxsmiles,
            atom_mapped=request.atom_map,
            rbi=balance_indices.rbi,
            pbi=balance_indices.pbi,
            tbi=balance_indices.tbi,
        )
    except Exception:
        logger.error("Error computing balance indices.", exc_info=True)
        raise HTTPException(status_code=400, detail="Error computing balance indices")


@reaction_utils_router.post("/atommap", summary="Map reactant atoms to product atoms using RXN4Chemistry RXNMapper")
async def get_atom_map(request: AtommapRequest) -> AtommapResponse:
    """
    Maps the given RXSMILES string using the RXNMapper and returns the mapped RXSMILES and confidence score.
    This function utilizes the [RXN4Chemistry RXNMapper package](https://github.com/rxn4chemistry/rxnmapper).

    Example input:
        { "rxsmiles": "CCO.CC(=O)O>>CC(=O)OCC.O" }
    """
    return get_atom_mapping_rxnmapper(request)


# TODO - Determine the purpose of this endpoint
@reaction_utils_router.get("/rxn_predict", summary="RXN prediction using TBD")
async def rxn_predict():
    """
    Predicts the reaction from RXNO.
    """
    # TODO - Copy /rxn_predict code from ASPIRE Prototype
    raise HTTPException(status_code=501, detail="Not implemented")


@reaction_utils_router.get("/rxn_cond_predict", summary="RXN condition prediction using TBD")
async def rxn_cond_predict():
    """
    Predicts the reaction condition from RXNO.
    """
    # TODO - Copy /rxn_predict code from ASPIRE Prototype
    raise HTTPException(status_code=501, detail="Not implemented")


@reaction_utils_router.post("/is_balanced", summary="Balanced reaction validation")
async def is_balanced(request: RxnIsBalancedRequest) -> RxnIsBalancedResponse:
    """
    Validates if a reaction is balanced.
    """
    try:
        is_balanced = is_reaction_balanced(request.rxsmiles)
        return RxnIsBalancedResponse(rxsmiles=request.rxsmiles, is_balanced=is_balanced)
    except Exception as e:
        logger.error(f"Error validating RXN balance: {str(e)}")
        raise HTTPException(status_code=500, detail="Error validating RXN balance")


@reaction_utils_router.post("/is_valid", summary="RXSMILES validation")
async def is_valid(requset: RxnIsValidRequest) -> RxnIsValidResponse:
    """
    Validates if a RXSMILES is valid.
    """
    try:
        is_valid = is_reaction_valid(requset.rxsmiles)
        return RxnIsValidResponse(rxsmiles=requset.rxsmiles, is_valid=is_valid)
    except Exception as e:
        logger.error(f"Error validating RXN: {str(e)}")
        raise HTTPException(status_code=500, detail="Error validating RXN")


@reaction_utils_router.get("/is_rxname_recognized", summary="RXName recognition")
async def is_rxname_recognized(rxname: str) -> bool:
    """
    Validates if a RXName is recognized.

    Valid examples: "Triflyloxy Menshutkin reaction", "Diazoalkane amination", "Formaldehyde reductive amination"
    Invalid examples: "Triflyloxy Menshutkin reaction 2", "Diazoalkane amination 2", "Invalid"
    """
    raise _hazelnut_not_available_http_exc()


@reaction_utils_router.get("/rxname_information", summary="Gets additional information about a RXName")
async def rxname_information(rxname: str):
    """
    Gets additional information about a RXName.

    Valid examples: "Triflyloxy Menshutkin reaction", "Diazoalkane amination", "Formaldehyde reductive amination"
    Invalid examples: "Triflyloxy Menshutkin reaction 2", "Diazoalkane amination 2", "Invalid"
    """
    raise _hazelnut_not_available_http_exc()


#########################
# NextMove Utils Router
#########################
nextmove_utils_router = APIRouter(prefix="/nextmove")


@nextmove_utils_router.post("/rxsmiles2rxno", summary="RXNO annotation with NameRXN from NextMove")
async def nextmove_rxno(request: HzRxnoRequest) -> HzRxnoResponse:
    """
    Retrieves the RXNO annotation of a reaction using NextMove's Hazelnut software.
    """
    raise _hazelnut_not_available_http_exc()


@nextmove_utils_router.post("/nm_atommap", summary="NextMove atom mapping")
async def nextmove_atommap(request: HzAtommapRequest) -> HzAtommapResponse:
    """
    Retrieves the atom map of a reaction using NextMove's Hazelnut software.
    """
    raise _hazelnut_not_available_http_exc()


@nextmove_utils_router.post("/normalize_role", summary="Normalize reaction role assignment")
async def nextmove_normalize_role(request: HzNormalizeRoleRequest) -> HzNormalizeRoleResponse:
    """
    Retrieves the normalized reaction role assignment using NextMove's Hazelnut software. Reaction must include atommapping.
    """
    raise _hazelnut_not_available_http_exc()


@nextmove_utils_router.post("/nm_map_and_normalize", summary="Provides atommapping and normalize role of reaction with NextMove")
async def nextmove_map_and_normalize(request: HzAtommapRequest) -> HzNormalizeRoleResponse:
    """
    Retrieves the atom map of a reaction then normalizes the roles using NextMove's Hazelnut software.
    """
    raise _hazelnut_not_available_http_exc()


#########################
# SVG Utils Router
#########################
svg_utils_router = APIRouter(prefix="/svg")


@svg_utils_router.post("/rxnsmiles2svg", summary="Convert RXN Smiles to SVG")
async def rxnsmiles2svg(request: Rxsmiles2SVGRequest) -> Rxsmiles2SVGResponse:
    """
    Converts RXN Smiles to SVG utilizing the RDKIT library.

    **Depiction modes:**
    - `simple`: Simple depiction, no atom mapping or highlighting.
    - `atom_map`: Atom map depiction.
    - `highlight_wo_indices`: Highlight without indices.
    - `highlight_with_indices`: Highlight with indices.

    **Args:**
    - `rxsmiles` (str): The RXN Smiles string to convert to SVG.
    - `depiction_mode` (str): The depiction mode to use. Options: `simple`, `atom_map`, `highlight_wo_indices`, `highlight_with_indices`.
    - `monochrome_atoms` (bool): If True, the atoms will be black.
    - `width` (int): The width of the SVG image.
    - `height` (int): The height of the SVG image.
    - `base64_encode` (bool): If True, the SVG will be base64 encoded.

    **Returns:**
    - `Rxsmiles2SVGResponse`: The RXN Smiles and the SVG image. Object with 'rxsmiles' and 'svg' keys.
    """
    try:
        svg = rxsmiles_to_svg(
            request.rxsmiles,
            depiction_mode=request.depiction_mode,
            monochrome_atoms=request.monochrome_atoms,
            img_width=request.width,
            img_height=request.height,
        )
        svg = svg.replace('"', "'")
        if request.base64_encode:
            svg = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
        return Rxsmiles2SVGResponse(rxsmiles=request.rxsmiles, svg=svg, base64_encoded=request.base64_encode)
    except Exception:
        logger.error("Error converting RXN to SVG:", exc_info=True)
        raise HTTPException(status_code=500, detail="Error converting RXN to SVG")


@svg_utils_router.post("/simple_rxnsmiles2svg", summary="Convert SMILES to SVG using ASKCOS method")
async def simple_rxnsmiles2svg(request: SimpleRxsmiles2SVGRequest) -> Rxsmiles2SVGResponse:
    """
    Converts RXN Smiles to SVG utilizing the ASKCOS method.

    **Args:**
    - `rxsmiles` (str): The RXN Smiles string to convert to SVG.
    - `base64_encode` (bool): If True, the SVG will be base64 encoded.
    - `retro` (bool): If True, the reaction is a retro-synthesis reaction.
    - `highlight_atoms` (bool): If True, the atoms will be highlighted in the reaction.
    - `show_atom_indices` (bool): If True, the atom indices will be shown in the reaction.

    **Returns:**
    - `Rxsmiles2SVGResponse`: The RXN Smiles and the SVG image. Object with 'rxsmiles' and 'svg' keys.
    """
    try:
        svg = reaction_smiles_to_image(request.rxsmiles, align=False, transparent=False, kekulize_mols=True,
                                       highlight=request.highlight_atoms, retro=request.retro, show_atom_indices=request.show_atom_indices)
        svg = svg.replace('"', "'")
    except Exception:
        svg = """
            <svg width="450" height="75" xmlns="http://www.w3.org/2000/svg">
            <rect width="100%" height="100%" fill="white" />
            <text x="10" y="50" font-size="32" fill="black">Unable to generate reaction SVG</text>
            </svg>
            """.strip()
        
    if request.base64_encode:
        svg = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
    return Rxsmiles2SVGResponse(rxsmiles=request.rxsmiles, svg=svg, base64_encoded=request.base64_encode)


# Include the routers in the main router
reaction_utils_router.include_router(nextmove_utils_router)
reaction_utils_router.include_router(svg_utils_router)
