from pydantic import BaseModel, Field


class HzRxnoRequest(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be mapped",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )


class HzRxnoResponse(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be mapped",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )
    rxname: str = Field(description="The name or identifier of the reaction.", examples=["2.6.8 O-Acetylatio"])
    rxclass: str = Field(description="The classification of the reaction, often indicating its type or category.", examples=["2.6 O-acylation to ester"])
    rxclass_rxno: str = Field(description="A unique identifier linking the reaction class with an RXNO ontology term.", examples=["MOP:0003479"])
    rxname_rxno: str = Field(description="A unique identifier linking the reaction name with an RXNO ontology term.", examples=[""])


class HzAtommapRequest(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be mapped",
        examples=["CCO.CC(=O)O>>CC(=O)OCC.O"],
    )
    maptype: str = Field(
        default="COMPLETER",
        description="NextMove Map type: 'MATCHING', 'COMPLETE', or 'COMPLETER'",
        examples=["COMPLETER"],
    )


class HzAtommapResponse(BaseModel):
    rxsmiles: str = Field(
        default="",
        title="RXSMILES",
        description="The RXSMILES string to be mapped",
        examples=[],
    )
    mapped_rxsmiles: str = Field(
        default="",
        title="Mapped RXSMILES",
        description="The mapped RXSMILES string",
        examples=[],
    )
    maptype: str = Field(
        description="NextMove Map type: 'MATCHING', 'COMPLETE', or 'COMPLETER'",
        examples=["COMPLETER"],
    )


class HzNormalizeRoleRequest(BaseModel):
    mapped_rxsmiles: str = Field(
        description="The reaction SMILES string with atom mappings to be normalized.",
        examples=["[CH3:1][CH2:2][OH:3].[CH3:4]C:5=[O:7]>>[CH3:1][CH2:2][O:3]C:5=[O:7].[OH2:6]"],
    )


class HzNormalizeRoleResponse(BaseModel):
    mapped_rxsmiles: str = Field(description="The original reaction SMILES string with atom mappings.", examples=[])
    role_normalized_rxsmiles: str = Field(description="The reaction SMILES string with normalized roles applied to atom mappings.", examples=[])
