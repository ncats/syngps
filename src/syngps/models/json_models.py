from enum import Enum
from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, Field, computed_field
from syngps.models import SynthGraph, SynthRoute, SynthRoutePreDefined
from syngps.utils import RxnSvgDepictionMode

##########################
# Reaction Models
##########################


class RxsmilesRequest(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be parsed",
        examples=[
            "[O:1]=[C:2]1[C:6]2([CH2:11][CH2:10][NH:9][CH2:8][CH2:7]2)[N:5]([C:12]2[CH:17]=[CH:16][CH:15]=[CH:14][CH:13]=2)[CH2:4][N:3]1[CH2:18][C:19]1[CH:31]=[CH:30][CH:29]=[CH:28][C:20]=1[C:21]([O:23][C:24]([CH3:27])([CH3:26])[CH3:25])=[O:22].[I-].[Na+].C(=O)([O-])[O-].[K+].[K+].Cl[CH2:41][CH2:42][CH2:43][N:44]1[C:52]2[C:47](=[CH:48][CH:49]=[CH:50][CH:51]=2)[C:46]([CH3:54])([CH3:53])[C:45]1=[O:55]>CC(=O)CC>[CH3:54][C:46]1([CH3:53])[C:47]2[C:52](=[CH:51][CH:50]=[CH:49][CH:48]=2)[N:44]([CH2:43][CH2:42][CH2:41][N:9]2[CH2:8][CH2:7][C:6]3([N:5]([C:12]4[CH:13]=[CH:14][CH:15]=[CH:16][CH:17]=4)[CH2:4][N:3]([CH2:18][C:19]4[CH:31]=[CH:30][CH:29]=[CH:28][C:20]=4[C:21]([O:23][C:24]([CH3:27])([CH3:25])[CH3:26])=[O:22])[C:2]3=[O:1])[CH2:11][CH2:10]2)[C:45]1=[O:55] |f:1.2,3.4.5|"
        ],
        min_length=1,
    )


class NormalizeRoleRequest(RxsmilesRequest):
    atom_map: bool = Field(
        default=False,
        title="Atom Map",
        description="Whether to add atom map to the RXSMILES before normalization",
        examples=[False],
    )


class NormalizeRoleResponse(RxsmilesRequest):
    original_rxsmiles: str = Field(
        default="",
        title="Original RXSMILES",
        description="The RXSMILES string to be parsed",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The normalized RXSMILES string",
        examples=[
            "[O:1]=[C:2]1[C:6]2([CH2:11][CH2:10][NH:9][CH2:8][CH2:7]2)[N:5]([C:12]2[CH:17]=[CH:16][CH:15]=[CH:14][CH:13]=2)[CH2:4][N:3]1[CH2:18][C:19]1[CH:31]=[CH:30][CH:29]=[CH:28][C:20]=1[C:21]([O:23][C:24]([CH3:27])([CH3:26])[CH3:25])=[O:22].Cl[CH2:41][CH2:42][CH2:43][N:44]1[C:52]2[C:47](=[CH:48][CH:49]=[CH:50][CH:51]=2)[C:46]([CH3:54])([CH3:53])[C:45]1=[O:55]>[K+].[K+].C(=O)([O-])[O-].CC(=O)CC.[Na+].[I-]>[CH3:54][C:46]1([CH3:53])[C:47]2[C:52](=[CH:51][CH:50]=[CH:49][CH:48]=2)[N:44]([CH2:43][CH2:42][CH2:41][N:9]2[CH2:8][CH2:7][C:6]3([N:5]([C:12]4[CH:13]=[CH:14][CH:15]=[CH:16][CH:17]=4)[CH2:4][N:3]([CH2:18][C:19]4[CH:31]=[CH:30][CH:29]=[CH:28][C:20]=4[C:21]([O:23][C:24]([CH3:27])([CH3:25])[CH3:26])=[O:22])[C:2]3=[O:1])[CH2:11][CH2:10]2)[C:45]1=[O:55] |f:2.3.4,6.7|"
        ],
    )
    atom_mapped: bool = Field(
        default=False,
        title="Atom Map",
        description="Whether to add atom map to the RXSMILES before normalization",
        examples=[False],
    )


class ComputeBalanceIndicesRequest(RxsmilesRequest):
    atom_map: bool = Field(
        default=False,
        title="Atom Map",
        description="Whether to add atom map to the RXSMILES before normalization",
        examples=[False],
    )


class ComputeBalanceIndicesResponse(RxsmilesRequest):
    rxsmiles: str = Field(default="", title="Original RXSMILES", description="The RXSMILES string to be parsed", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])
    atom_mapped: bool = Field(
        default=False,
        title="Atom Map",
        description="Whether to add atom map to the RXSMILES before computation",
        examples=[False],
    )
    rbi: float = Field(default=0.0, title="Reactant Balance Index", description="The reactant balance index of the RXSMILES", examples=[0.0])
    pbi: float = Field(default=0.0, title="Product Balance Index", description="The product balance index of the RXSMILES", examples=[0.0])
    tbi: float = Field(default=0.0, title="Total Balance Index", description="The total balance index of the RXSMILES", examples=[0.0])


class ParseRxsmilesResponse(BaseModel):
    rxsmiles: str = Field(default="", title="RXSMILES", description="The original RXSMILES that was parsed", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])
    fragment_groups: Optional[List[List[int]]] = Field(
        default=None, title="Fragment Groups", description="The fragment groups in the RXSMILES", examples=[None]
    )
    reactants: List[str] = Field(default=[], title="Reactants", description="The reactants in the reaction", examples=["CCO", "CC(=O)O"])
    reagents: List[str] = Field(default=[], title="Reagents", description="The reagents in the reaction", examples=[])
    products: List[str] = Field(default=[], title="Products", description="The products in the reaction", examples=["CC(=O)OCC", "O"])


class AtommapRequest(BaseModel):
    smiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be mapped",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )


class AtommapResponse(BaseModel):
    original_rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be mapped",
        examples=["[C:4]([CH3:5])(=[O:6])[OH:7].[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][O:3][C:4]([CH3:5])=[O:6].[OH2:7]"],
    )
    mapped_rxsmiles: str = Field(
        default="",
        title="Mapped RXSMILES",
        description="The mapped RXSMILES string",
        examples=["[C:4]([CH3:5])(=[O:6])[OH:7].[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][O:3][C:4]([CH3:5])=[O:6].[OH2:7]"],
    )
    confidence: float = Field(
        default=0.0,
        title="Confidence",
        description="The confidence score of the mapping",
        examples=[0.8266466076655647],
    )


class RxnIsBalancedRequest(BaseModel):
    rxsmiles: str = Field(default="", title="RXN", description="The RXN string to be validated", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])


class RxnIsBalancedResponse(BaseModel):
    rxsmiles: str = Field(default="", title="RXN", description="The RXN string to be validated", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])
    is_balanced: bool = Field(default=False, title="Is Balanced", description="Whether the RXN is balanced", examples=[True])


class RxnIsValidRequest(BaseModel):
    rxsmiles: str = Field(default="", title="RXN", description="The RXN string to be validated", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])


class RxnIsValidResponse(BaseModel):
    rxsmiles: str = Field(default="", title="RXN", description="The RXN string to be validated", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])
    is_valid: bool = Field(default=False, title="Is Valid", description="Whether the RXN is valid", examples=[True])


class Rxsmiles2SVGRequest(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXN",
        description="The RXN string to be converted to SVG",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )
    depiction_mode: RxnSvgDepictionMode = Field(
        default=RxnSvgDepictionMode.SIMPLE.value,
        title="RXN Depiction Mode",
        description="The depiction mode of the RXN. Can be simple, atom_map, highlight_wo_indices or highlight_with_indices",
        examples=["simple", "atom_map", "highlight_wo_indices", "highlight_with_indices"],
    )
    monochrome_atoms: bool = Field(
        default=False,
        title="Monochrome Atoms",
        description="Whether to render the atoms in black or color",
        examples=[False],
    )
    width: int = Field(
        default=1200,
        title="Width",
        description="The width of the SVG",
        examples=[1200],
    )
    height: int = Field(
        default=400,
        title="Height",
        description="The height of the SVG",
        examples=[400],
    )
    base64_encode: bool = Field(
        default=False,
        title="Base64 Encode",
        description="Whether to base64 encode the SVG",
        examples=[False],
    )


class SimpleRxsmiles2SVGRequest(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXN",
        description="The RXN string to be converted to SVG",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )
    retro: bool = Field(
        default=False,
        title="Retro Depiction",
        description="Whether this is a retrosynthetic reaction",
        examples=[False],
    )
    highlight_atoms: bool = Field(
        default=True,
        title="Highlight",
        description="Whether to atom highlight atoms the reaction",
        examples=[True],
    )
    base64_encode: bool = Field(
        default=False,
        title="Base64 Encode",
        description="Whether to base64 encode the SVG",
        examples=[False],
    )


class Rxsmiles2SVGResponse(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXN",
        description="The RXN string to be converted to SVG",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )
    svg: str = Field(
        default="",
        title="SVG",
        description="The SVG representation of the RXN",
        examples=['<svg xmlns="http://www.w3.org/2000/svg">EXAMPLE</svg>'],
    )
    base64_encoded: bool = Field(
        default=False,
        title="Base64 Encoded",
        description="Whether the SVG is base64 encoded",
        examples=[False],
    )


class Rxid2SVGRequest(BaseModel):
    rxid: str = Field(
        default="",
        title="RXID",
        description="The RXID string to be converted to SVG",
        examples=["ASPIRE-00017180000384"],
    )
    depiction_mode: RxnSvgDepictionMode = Field(
        default=RxnSvgDepictionMode.SIMPLE,
        title="RXN Depiction Mode",
        description="The depiction mode of the RXN. Can be simple, atom_map, highlight_wo_indices or highlight_with_indices",
        examples=["simple", "atom_map", "highlight_wo_indices", "highlight_with_indices"],
    )
    monochrome_atoms: bool = Field(
        default=False,
        title="Monochrome Atoms",
        description="Whether to render the atoms in black or color",
        examples=[False],
    )
    width: int = Field(
        default=1200,
        title="Width",
        description="The width of the SVG",
        examples=[1200],
    )
    height: int = Field(
        default=400,
        title="Height",
        description="The height of the SVG",
        examples=[400],
    )
    base64_encode: bool = Field(
        default=False,
        title="Base64 Encode",
        description="Whether to base64 encode the SVG",
        examples=[False],
    )


class Rxid2SVGResponse(BaseModel):
    rxid: str = Field(
        default="",
        title="RXID",
        description="The RXID string to be converted to SVG",
        examples=["ASPIRE-00017180000384"],
    )
    svg: str = Field(
        default="",
        title="SVG",
        description="The SVG representation of the RXN",
        examples=['<svg xmlns="http://www.w3.org/2000/svg">EXAMPLE</svg>'],
    )
    base64_encoded: bool = Field(
        default=False,
        title="Base64 Encoded",
        description="Whether the SVG is base64 encoded",
        examples=[False],
    )


class RxnClassifyRequest(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be classified",
        examples=["Cc1ccc2[nH]cc(CCN)c2c1.O=CC1CC1>>Cc1ccc2[nH]cc(CCNCC3CC3)c2c1"],
    )
    return_probabilities: bool = Field(
        default=True,
        title="Return Probabilities",
        description="Whether to return class probabilities",
        examples=[True],
    )
    return_top_n: int = Field(
        default=5,
        title="Return Top N",
        description="The number of top predictions to return",
        examples=[5],
    )


class RxnClassifyResponse(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be classified",
        examples=["Cc1ccc2[nH]cc(CCN)c2c1.O=CC1CC1>>Cc1ccc2[nH]cc(CCNCC3CC3)c2c1"],
    )
    classes: List[str] = Field(
        default=[],
        title="Classes",
        description="The predicted classes of the RXN",
        examples=[
            [
                "1.2.1 Aldehyde reductive amination",
                "1.2.5 Ketone reductive amination",
                "0.0.0 Other reaction",
                "5.1.1 N-Boc protection",
                "4.2.1 1,2,4-Oxadiazol-5-one synthesis",
            ]
        ],
    )
    probabilities: List[float] | None = Field(
        default=None,
        title="Probabilities",
        description="The predicted probabilities of the RXN",
        examples=[[0.9958903958557691, 0.004059463689513187, 0.000011208619127173303, 0.000007613874398914294, 0.0000041039589689487136]],
    )


##########################
# Substance Models
##########################


class SubstanceSmiles2InchikeyRequest(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to InChI key",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )


class SubstanceSmiles2InchikeyResponse(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to InChI key",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )
    inchikey: str = Field(
        default="",
        title="InChI Key",
        description="The InChI key of the substance",
        examples=["AEASCKUEJBJPEL-UHFFFAOYSA-N"],
    )


class SubstanceSmiles2BMSRequest(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to BMS",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )


class SubstanceSmiles2BMSResponse(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to BMS",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )
    bms: str = Field(
        default="",
        title="BMS",
        description="The Bond Manipulation System of the substance",
        examples=["O=C1CCC(=O)C1c1ccccc1"],
    )


class SubstanceSmiles2SDFRequest(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to SDF",
        examples=["CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"],
    )
    toKekulize: Optional[bool] = Field(
        default=False,
        title="Kekulization",
        description="Flag for whether the molecule should be kekulized",
        examples=[False],
    )


class SubstanceSmiles2SDFResponse(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to SDF",
        examples=["CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"],
    )
    sdf: str = Field(
        default="",
        title="SDF",
        description="The Structural Data File of the substance",
        examples=[
            "\n"
            + "  Mrv2301 09232414262D          \n"
            + "\n"
            + " 15 15  0  0  1  0            999 V2000\n"
            + "    0.0000    3.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    2.4750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -1.2375    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -2.4750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145   -1.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  1  0  0  0  0\n"
            + "  2  3  1  0  0  0  0\n"
            + "  2  4  1  0  0  0  0\n"
            + "  4  5  1  0  0  0  0\n"
            + "  5  6  4  0  0  0  0\n"
            + "  6  7  4  0  0  0  0\n"
            + "  7  8  4  0  0  0  0\n"
            + "  8  9  4  0  0  0  0\n"
            + "  9 10  4  0  0  0  0\n"
            + "  5 10  4  0  0  0  0\n"
            + " 11  8  1  0  0  0  0\n"
            + " 11 12  1  1  0  0  0\n"
            + " 11 13  1  0  0  0  0\n"
            + " 13 14  2  0  0  0  0\n"
            + " 13 15  1  0  0  0  0\n"
            + "M  END\n"
            + "$$$$"
        ],
    )


class SubstanceSDF2SmilesRequest(BaseModel):
    sdf: str = Field(
        default="",
        title="SDF",
        description="The Structural Data File to be converted to SMILES string",
        examples=[
            "\n"
            + "  Mrv2301 09232414262D          \n"
            + "\n"
            + " 15 15  0  0  1  0            999 V2000\n"
            + "    0.0000    3.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    2.4750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -1.2375    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -2.4750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145   -1.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  1  0  0  0  0\n"
            + "  2  3  1  0  0  0  0\n"
            + "  2  4  1  0  0  0  0\n"
            + "  4  5  1  0  0  0  0\n"
            + "  5  6  4  0  0  0  0\n"
            + "  6  7  4  0  0  0  0\n"
            + "  7  8  4  0  0  0  0\n"
            + "  8  9  4  0  0  0  0\n"
            + "  9 10  4  0  0  0  0\n"
            + "  5 10  4  0  0  0  0\n"
            + " 11  8  1  0  0  0  0\n"
            + " 11 12  1  1  0  0  0\n"
            + " 11 13  1  0  0  0  0\n"
            + " 13 14  2  0  0  0  0\n"
            + " 13 15  1  0  0  0  0\n"
            + "M  END\n"
            + "$$$$"
        ],
    )


class SubstanceSDF2SmilesResponse(BaseModel):
    sdf: str = Field(
        default="",
        title="SDF",
        description="The Structural Data File to be converted to SMILES string",
        examples=[
            "\n"
            + "  Mrv2301 09232414262D          \n"
            + "\n"
            + " 15 15  0  0  1  0            999 V2000\n"
            + "    0.0000    3.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    2.4750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -1.2375    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -2.4750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145   -1.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  1  0  0  0  0\n"
            + "  2  3  1  0  0  0  0\n"
            + "  2  4  1  0  0  0  0\n"
            + "  4  5  1  0  0  0  0\n"
            + "  5  6  4  0  0  0  0\n"
            + "  6  7  4  0  0  0  0\n"
            + "  7  8  4  0  0  0  0\n"
            + "  8  9  4  0  0  0  0\n"
            + "  9 10  4  0  0  0  0\n"
            + "  5 10  4  0  0  0  0\n"
            + " 11  8  1  0  0  0  0\n"
            + " 11 12  1  1  0  0  0\n"
            + " 11 13  1  0  0  0  0\n"
            + " 13 14  2  0  0  0  0\n"
            + " 13 15  1  0  0  0  0\n"
            + "M  END\n"
            + "$$$$"
        ],
    )
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string of the substance",
        examples=["CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"],
    )


class SubstanceSDFIsValidRequest(BaseModel):
    sdf: str = Field(
        default="",
        title="SDF",
        description="The Structural Data File to be checked for validity against RDKit",
        examples=[
            "\n"
            + "  Mrv2301 09232414262D          \n"
            + "\n"
            + " 15 15  0  0  1  0            999 V2000\n"
            + "    0.0000    3.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    2.4750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -1.2375    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -2.4750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145   -1.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  1  0  0  0  0\n"
            + "  2  3  1  0  0  0  0\n"
            + "  2  4  1  0  0  0  0\n"
            + "  4  5  1  0  0  0  0\n"
            + "  5  6  4  0  0  0  0\n"
            + "  6  7  4  0  0  0  0\n"
            + "  7  8  4  0  0  0  0\n"
            + "  8  9  4  0  0  0  0\n"
            + "  9 10  4  0  0  0  0\n"
            + "  5 10  4  0  0  0  0\n"
            + " 11  8  1  0  0  0  0\n"
            + " 11 12  1  1  0  0  0\n"
            + " 11 13  1  0  0  0  0\n"
            + " 13 14  2  0  0  0  0\n"
            + " 13 15  1  0  0  0  0\n"
            + "M  END\n"
            + "$$$$"
        ],
    )


class SubstanceSDFIsValidResponse(BaseModel):
    sdf: str = Field(
        default="",
        title="SDF",
        description="The Structural Data File to be checked for validity against RDKit",
        examples=[
            "\n"
            + "  Mrv2301 09232414262D          \n"
            + "\n"
            + " 15 15  0  0  1  0            999 V2000\n"
            + "    0.0000    3.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    2.4750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.7145   -1.2375    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"
            + "    1.4289   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.0000   -2.4750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "   -0.7145   -1.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  1  0  0  0  0\n"
            + "  2  3  1  0  0  0  0\n"
            + "  2  4  1  0  0  0  0\n"
            + "  4  5  1  0  0  0  0\n"
            + "  5  6  4  0  0  0  0\n"
            + "  6  7  4  0  0  0  0\n"
            + "  7  8  4  0  0  0  0\n"
            + "  8  9  4  0  0  0  0\n"
            + "  9 10  4  0  0  0  0\n"
            + "  5 10  4  0  0  0  0\n"
            + " 11  8  1  0  0  0  0\n"
            + " 11 12  1  1  0  0  0\n"
            + " 11 13  1  0  0  0  0\n"
            + " 13 14  2  0  0  0  0\n"
            + " 13 15  1  0  0  0  0\n"
            + "M  END\n"
            + "$$$$"
        ],
    )
    valid: bool = Field(
        default=True,
        title="Validity",
        description="Whether the SDF is parseable by RDKit",
        examples=[True],
    )


class SubstanceSmiles2SVGRequest(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to SVG",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )
    monochrome_atoms: bool = Field(
        default=False,
        title="Monochrome Atoms",
        description="Whether to render the atoms in black or color",
        examples=[False],
    )
    # kekulize_mol: bool = Field(
    #     default=False,
    #     title="Kekulize Molecule",
    #     description="Whether to kekulize the molecule",
    #     examples=[False],
    # )
    width: int = Field(
        default=300,
        title="Width",
        description="The width of the SVG",
        examples=[300],
    )
    height: int = Field(
        default=300,
        title="Height",
        description="The height of the SVG",
        examples=[300],
    )
    base64_encode: bool = Field(
        default=False,
        title="Base64 Encode",
        description="Whether to base64 encode the SVG",
        examples=[False],
    )


class SimpleSubstanceSmiles2SVGRequest(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to SVG",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )
    base64_encode: bool = Field(
        default=False,
        title="Base64 Encode",
        description="Whether to base64 encode the SVG",
        examples=[False],
    )


class SubstanceSmiles2SVGResponse(BaseModel):
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string to be converted to SVG",
        examples=["Cc1cc(Br)cc(C)c1C1C(=O)CCC1=O"],
    )
    svg: str = Field(
        default="",
        title="SVG",
        description="The SVG representation of the substance",
        examples=['<svg xmlns="http://www.w3.org/2000/svg">EXAMPLE</svg>'],
    )
    base64_encoded: bool = Field(
        default=False,
        title="Base64 Encoded",
        description="Whether the SVG is base64 encoded",
        examples=[False],
    )


class SubstanceInchikey2SVGRequest(BaseModel):
    inchikey: str = Field(
        default="",
        title="InChI Key",
        description="The InChI key to be converted to SVG",
        examples=["AEASCKUEJBJPEL-UHFFFAOYSA-N"],
    )
    monochrome_atoms: bool = Field(
        default=False,
        title="Monochrome Atoms",
        description="Whether to render the atoms in black or color",
        examples=[False],
    )
    # kekulize_mol: bool = Field(
    #     default=False,
    #     title="Kekulize Molecule",
    #     description="Whether to kekulize the molecule",
    #     examples=[False],
    # )
    width: int = Field(
        default=300,
        title="Width",
        description="The width of the SVG",
        examples=[300],
    )
    height: int = Field(
        default=300,
        title="Height",
        description="The height of the SVG",
        examples=[300],
    )
    base64_encode: bool = Field(
        default=False,
        title="Base64 Encode",
        description="Whether to base64 encode the SVG",
        examples=[False],
    )


class SubstanceInchikey2SVGResponse(BaseModel):
    inchikey: str = Field(
        default="",
        title="InChI Key",
        description="The InChI key to be converted to SVG",
        examples=["AAAJPELFCSVKMB-HWKANZROSA-N"],
    )
    smiles: str = Field(
        default="",
        title="SMILES",
        description="The SMILES string of the substance",
        examples=["CN1CCN(CC/C=C/c2cc3ncc(C#N)c(Nc4ccc(Sc5nccn5C)c(Cl)c4)c3s2)CC1"],
    )
    svg: str = Field(
        default="",
        title="SVG",
        description="The SVG representation of the substance",
        examples=['<svg xmlns="http://www.w3.org/2000/svg">EXAMPLE</svg>'],
    )
    base64_encoded: bool = Field(
        default=False,
        title="Base64 Encoded",
        description="Whether the SVG is base64 encoded",
        examples=[False],
    )


############################
# ASPIRE Inventory Models
############################


# Define Inventory Status enum
class InventoryStatus(str, Enum):
    IN_STOCK = "In Stock - Stereo Match"  # Exact Inchikey match
    IN_STOCK_NS = "In Stock - Non-Stereo Match"  # NS Inchikey match
    NOT_AVAILABLE = "Not Available"
    UNKNOWN = "Unknown"


class InventoryStatusRequest(BaseModel):
    inchikeys: List[str] = Field(
        default=[],
        title="Inchikeys to retrieve inventory statuses for",
        description="List of InChI keys to retrieve inventory statuses for",
        examples=[["CKMVPXUSFQBAIJ-UHFFFAOYSA-N", "AAAJPELFCSVKMB-HWKANZROSA-N"]],
    )


class InventoryStatusResponse(BaseModel):
    inchikeys: List[str] = Field(
        default=[],
        title="Inchikeys to retrieve inventory statuses for",
        description="List of InChI keys to retrieve inventory statuses for",
        examples=[["CKMVPXUSFQBAIJ-UHFFFAOYSA-N", "AAAJPELFCSVKMB-HWKANZROSA-N"]],
    )
    inventory_statuses: Dict[str, InventoryStatus] = Field(
        default={},
        title="Inventory Statuses",
        description="Dictionary of InChI keys to inventory statuses",
        examples=[{"CKMVPXUSFQBAIJ-UHFFFAOYSA-N": InventoryStatus.IN_STOCK, "AAAJPELFCSVKMB-HWKANZROSA-N": InventoryStatus.NOT_AVAILABLE}],
    )


############################
# ASPIRE Synthesis Models
############################


class EvidenceBasedRouteDetails(BaseModel):
    query_type: str = Field(default="shortest_path", examples=["shortest_path"])
    top_n_routes: int = Field(default=3, examples=[3], ge=1)


class PredictedRouteDetails(BaseModel):
    source: Literal["ASKCOS v2"]
    max_routes: int = Field(default=3, examples=[3])


class SynthesisRoutesRequest(BaseModel):
    target_molecule_inchikey: Optional[str] = Field(examples=["BACODVKOMBXPLL-UHFFFAOYSA-N"], default=None)
    target_molecule_smiles: Optional[str] = Field(examples=["O=C(O)c1cc(-c2ccccc2)cc2cc[nH]c12"], default=None)
    reaction_steps: int = Field(default=2, examples=[2], gt=0)
    include_svgs: bool = Field(default=False, examples=[False])
    include_evidence_routes: bool = Field(default=True, examples=[True])
    evidence_options: EvidenceBasedRouteDetails = Field()
    include_evidence_synth_graph: bool = Field(default=True, examples=[True])
    include_predicted_routes: bool = Field(default=True, examples=[True])
    prediction_options: PredictedRouteDetails = Field()
    include_predictive_synth_graph: bool = Field(default=True, examples=[True])


class SynthesisRoutesResponse(BaseModel):
    target_molecule_inchikey: str = Field(examples=["BACODVKOMBXPLL-UHFFFAOYSA-N"], min_length=1)
    target_molecule_smiles: str = Field(examples=["O=C(O)c1cc(-c2ccccc2)cc2cc[nH]c12"], min_length=1)
    reaction_steps: int = Field(default=2, examples=[2], gt=0)
    evidence_routes_success: bool = Field(default=False, examples=[True])
    predicted_routes_success: bool = Field(default=False, examples=[True])
    routes: List[SynthRoute | SynthRoutePreDefined] = Field(default=[])
    evidence_synth_graph: Optional[SynthGraph] = Field(default=None)
    predictive_synth_graph: Optional[SynthGraph] = Field(default=None)

    @computed_field
    def num_routes(self) -> int:
        return len(self.routes)
    
    @computed_field
    def num_evidence_routes(self) -> int:
        return len([route for route in self.routes if route.predicted == False])
    
    @computed_field
    def num_predicted_routes(self) -> int:
        return len([route for route in self.routes if route.predicted == True])
