import math
from pstats import Stats
from typing import Any, Dict, List, Optional, Union

from networkx import DiGraph
from pydantic import BaseModel, ConfigDict, Field, computed_field, model_validator
from utils import (
    graph_to_json,
    graph_to_json_edges,
    graph_to_json_nodes,
)


class ProfiledResult(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    profiling_stats: Optional[Stats] = Field(default=None, exclude=True)

    class CustomConfig:
        @staticmethod
        def schema_extra(schema: Dict[str, Any], model: BaseModel) -> None:
            if "properties" in schema:
                schema["properties"].pop("profiling_stats", None)


class TimeInfo(BaseModel):
    synth_graph_time: Optional[float] = Field(default=None)
    enumerate_routes_time: Optional[float] = Field(default=None)
    aggregated_yield_time: Optional[float] = Field(default=None)
    rank_routes_time: Optional[float] = Field(default=None)


class TopologyInfo(BaseModel):
    num_edges: Optional[int] = Field(default=None)
    num_nodes: Optional[int] = Field(default=None)
    num_routes: Optional[int] = Field(default=None)
    num_combination_graphs: Optional[int] = Field(default=None)
    num_route_candidates: Optional[int] = Field(default=None)


class Location(BaseModel):
    smiles: str = Field(description="SMILES representation of the substance")
    room: Optional[str] = Field(default=None, description="Room where the substance is located")
    position: Optional[str] = Field(default=None, description="Position in the room")
    quantity_weight: Optional[str] = Field(default=None, description="Quantity weight of the substance")
    unit: Optional[str] = Field(default=None, description="Unit of the quantity weight")


class Vendor(BaseModel):
    smiles: str = Field(description="SMILES representation of the substance")
    source: Optional[str] = Field(default=None, description="Source of the substance")
    ppg: Optional[str] = Field(default=None, description="Price per gram")
    lead_time: Optional[str] = Field(default=None, description="Lead time for availability")
    url: Optional[str] = Field(default=None, description="URL for the vendor")


class Inventory(BaseModel):
    available: bool = Field(description="Availability status in inventory")
    locations: List[Location] = Field(default_factory=list, description="List of locations where the substance is stored")


class CommercialAvailability(BaseModel):
    available: bool = Field(description="Availability status commercially")
    vendors: List[Vendor] = Field(default_factory=list, description="List of vendors offering the substance")


class Availability(BaseModel):
    inchikey: str = Field(description="InChIKey of the substance")
    inventory: Inventory = Field(description="Inventory availability details")
    commercial_availability: CommercialAvailability = Field(description="Commercial availability details")


class SynthGraphSearch(BaseModel):
    target_molecule_inchikey: str = Field(examples=["YTIQRXMAVJLXHH-UHFFFAOYSA-N"], min_length=1)
    reaction_steps: int = Field(default=2, examples=[2], description="Number of reaction steps to search away from target molecule", gt=0)
    query_type: str = Field(
        default="shortest_path",
        examples=["shortest_path"],
        description='Type of query to be used for extracting graph data, options are "shortest_path". "shortest_path" (the default) will use shortest path to truncate much of the graph searching to enable faster query time.',
    )
    leaves_as_sm: bool = Field(default=True, examples=[True], description="Whether to consider leaf substance nodes as starting materials")
    include_availability_info: bool = Field(default=False, examples=[False], description="Whether to include availability information in the response")
    annotate_reactions: bool = Field(default=False, examples=[False], description="Whether to annotate reactions in the graph with Hazelnut")
    graph_backend: str = Field(default="memgraph", examples=[
                               "memgraph", "neo4j"], description="Whether to use memgraph or neo4j backend. Neo4j backend may not be enabled for all deployments.")

    @property
    def search_depth(self) -> int:
        # Automatically calculate search_depth as 2 * reaction_steps
        return 2 * self.reaction_steps


class TopNYieldSearch(SynthGraphSearch):
    top_n_routes: int = Field(default=3, examples=[3], ge=1)
    include_route_candidates: bool = Field(default=False, examples=[False], description="Whether to return route candidates or not")
    include_combination_graphs: bool = Field(default=False, examples=[False], description="Whether to return combination graphs or not")
    synthesis_graph_json: Optional[Dict[str, Any]] = Field(default=None, description="Optional synthesis graph input", examples=[None])


class SynthGraphJsonInput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    nodes: List[Dict] = Field(default=None)
    edges: List[Dict] = Field(default=None)


class SynthRoute(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    route_candidate: DiGraph = Field(exclude=True)
    aggregated_yield: Optional[float] = Field(default=None)
    predicted: bool = Field(default=False)
    source: Optional[str] = Field(default=None)
    route_index: Optional[int] = Field(default=None)
    route_status: Optional[str] = Field(default=None)
    method: Optional[str] = Field(default=None)

    @computed_field
    def route_node_labels(self) -> List[str]:
        """Extract and return the node_label of all nodes in the route_candidate graph."""
        if not hasattr(self, "route_candidate") or not self.route_candidate:
            return []
        try:
            return [
                node_data.get("node_label")
                for _, node_data in self.route_candidate.nodes(data=True)
                if "node_label" in node_data
            ]
        except Exception:
            return []


class CombinationGraph(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    route_index: int = Field(default=None)
    route_candidate: DiGraph = Field(exclude=True)

    @computed_field
    def route_candidate_json(self) -> Dict[str, Any]:
        return graph_to_json(self.route_candidate)


class SynthRoutePreDefined(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    route_candidate: None = None
    aggregated_yield: Optional[float] = Field(default=None)
    predicted: bool = Field(default=False)
    source: Optional[str] = Field(default=None)
    route_index: Optional[int] = Field(default=None)
    route_status: Optional[str] = Field(default=None)
    method: Optional[str] = Field(default=None)
    route_node_labels: List[str] = Field(default=None)


class SynthView(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    view: DiGraph = Field(exclude=True)

    @computed_field
    def route_node_labels(self) -> List[str]:
        """Extract and return the node_label of all nodes in the route_candidate graph."""
        return [node_data.get("node_label") for _, node_data in self.view.nodes(data=True)]


class ProjectedRoute(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    route_index: int = Field(default=None)
    route: DiGraph = Field(exclude=True)
    aggregated_yield: Optional[float] = Field(default=None)

    @computed_field
    def route_json(self) -> dict:
        return graph_to_json(self.route)


class Substance(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    inchikey: str = Field()
    nsinchikey: Optional[str] = Field(default=None)
    canonical_smiles: str = Field(default=None)
    graph_node_id: Optional[str | int] = Field(default=None)
    srole: Optional[str] = Field(default=None)
    is_in_inventory: Optional[bool] = Field(default=None)
    is_in_uspto_full: bool = Field(default=False)
    is_in_savi_130k: bool = Field(default=False)
    is_in_aicp: bool = Field(default=False)
    product_in: Optional[List[str]] = Field(default=None)
    reactant_in: Optional[List[str]] = Field(default=None)
    reagent_in: Optional[List[str]] = Field(default=None)
    is_predicted: Optional[bool] = Field(default=None)
    is_evidence: Optional[bool] = Field(default=None)


class Reaction(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    # Define model fields
    rxid: str = Field()
    rxsmiles: str = Field()
    graph_node_id: Optional[str | int] = Field(default=None)
    eln_exp_name: Optional[str] = Field(default=None)
    eln_exp_id: Optional[str] = Field(default=None)
    eln_rxn_id: Optional[str] = Field(default=None)
    is_balanced: bool = Field(default=None)
    is_rxname_recognized: bool = Field(default=None)
    is_rxsmiles_valid: bool = Field(default=None)
    is_in_uspto_full: Optional[bool] = Field(default=None)
    is_in_savi_130k: Optional[bool] = Field(default=None)
    is_in_aicp: Optional[bool] = Field(default=None)
    rbi: Optional[float] = Field(default=None)
    pbi: Optional[float] = Field(default=None)
    tbi: Optional[float] = Field(default=None)
    yield_predicted: Optional[float] = Field(default=None)
    additional_data: Optional[Dict[str, Any]] = Field(default=None)
    products: Optional[List[Substance]] = Field(default=None)
    reactants: Optional[List[Substance]] = Field(default=None)
    reagents: Optional[List[Substance]] = Field(default=None)
    rxsmiles_original: Optional[str] = Field(default=None)
    patent_number: Optional[str] = Field(default=None, alias="PatentNumber", serialization_alias="patent_number")
    paragraph_number: Optional[Union[int, str]] = Field(default=None, alias="ParagraphNum", serialization_alias="paragraph_number")
    is_predicted: Optional[bool] = Field(default=None)
    is_evidence: Optional[bool] = Field(default=None)
    yield_empirical: Optional[float] = Field(default=None)

    @model_validator(mode="before")
    def replace_all_nans(cls, data):
        for key, value in data.items():
            if isinstance(value, float) and math.isnan(value):
                data[key] = None
        return data


class ReactionSource(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    # Define model fields
    rxid: str = Field()
    rxsmiles: str = Field()
    eln_exp_name: Optional[str] = Field(default=None)
    eln_exp_id: Optional[str] = Field(default=None)
    eln_rxn_id: Optional[str] = Field(default=None)
    is_in_uspto_full: Optional[bool] = Field(default=None)
    is_in_savi_130k: Optional[bool] = Field(default=None)
    is_in_aicp: Optional[bool] = Field(default=None)
    rxsmiles_original: Optional[str] = Field(default=None)
    patent_number: Optional[str] = Field(default=None, alias="PatentNumber", serialization_alias="patent_number")
    paragraph_number: Optional[Union[int, str]] = Field(default=None, alias="ParagraphNum", serialization_alias="paragraph_number")


class SynthGraph(ProfiledResult):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    synthesis_graph: DiGraph = Field(exclude=True)
    search_params: SynthGraphSearch | None = Field(default=None)
    target_molecule_node_id: str | int = Field(exclude=True)
    time_info: Optional[TimeInfo] = Field(default=None)
    summary: Optional[TopologyInfo] = Field(default=None)
    availability: Optional[List[Availability]] = Field(default=None, description="Availability information for substances")

    @computed_field
    def nodes(self) -> list:
        return graph_to_json_nodes(self.synthesis_graph)

    @computed_field
    def edges(self) -> list:
        return graph_to_json_edges(self.synthesis_graph)

    def hash_graph(self) -> str:
        return hash_graph(self.synthesis_graph)


class SynthRouteContainer(BaseModel):
    """Encapsulates a list of SynthRoute objects along with the number of subgraphs."""

    subgraphs: Optional[List[SynthRoute]] = Field(default_factory=list)

    @computed_field
    def num_subgraphs(self) -> int:
        """Returns the number of subgraphs, handling empty or None cases."""
        return len(self.subgraphs) if self.subgraphs else 0


class RouteAssemblyType(BaseModel):
    is_predicted: bool
    is_evidence: bool


class Node(BaseModel):
    node_label: str
    node_type: str
    uuid: str


class Provenance(BaseModel):
    is_in_aicp: bool
    is_in_uspto_full: bool
    is_in_savi_130k: bool
    patents: Optional[List[str]] = Field(default=None, description="List of patent strings")
    patent_paragraph_nums: Optional[List[int]] = Field(default=None, description="List of paragraph numbers corresponding to patents")


class SubstanceNode(Node):
    node_type: str = "substance"
    inchikey: str
    canonical_smiles: str
    srole: str
    route_assembly_type: RouteAssemblyType
    provenance: Provenance = Field(description="Provenance information for the substance node")


class Validation(BaseModel):
    is_balanced: bool = Field(default=None, description="Whether the reaction is balanced")
    is_rxname_recognized: Optional[bool] = Field(default=None, description="Whether the reaction name is recognized")
    is_valid: Optional[bool] = Field(default=True, description="Whether the validation is successful")


class YieldInfo(BaseModel):
    yield_predicted: float
    yield_score: int


class EvidenceProtocol(BaseModel):
    protocol_text: Dict[str, str] = Field(description="Map of keys (e.g., patents) to detailed protocol texts for the reaction evidence")


class ReactionNode(Node):
    node_type: str = Field(default="reaction", description="Type of the node, e.g., reaction")
    rxid: str = Field(default=None, description="Reaction ID")
    rxsmiles: str = Field(default=None, description="Reaction SMILES string")
    rxclass: str = Field(default=None, description="Reaction class")
    rxname: str = Field(default=None, description="Reaction name")
    original_rxsmiles: str = Field(default=None, description="Original reaction SMILES")
    yield_info: YieldInfo
    provenance: Provenance  # Contains is_in_aicp
    validation: Validation  # Contains is_balanced
    route_assembly_type: RouteAssemblyType
    evidence_protocol: Optional[Dict[str, EvidenceProtocol]] = Field(default=None, description="Evidence protocol for the reaction")
    evidence_conditions_info: Optional[Dict[str, Dict[str, str]]] = Field(
        default=None, description="Evidence conditions info for the reaction, mapping keys (e.g., patents) to dictionaries of condition key-value pairs"
    )
    predicted_conditions_info: Optional[Dict[str, Dict[str, Dict[str, str]]]] = Field(
        default=None,
        description="Predicted conditions info for the reaction, mapping methods to dictionaries of predictions, each containing condition key-value pairs",
    )


class Edge(BaseModel):
    start_node: str
    end_node: str
    edge_label: str
    edge_type: str
    provenance: Provenance  # Contains is_in_aicp
    uuid: str
    inchikey: str
    rxid: str
    route_assembly_type: RouteAssemblyType  # Using RouteAssemblyType instead of a dict


class SynthGraphJson(BaseModel):
    nodes: List[Union[SubstanceNode, ReactionNode]]
    edges: List[Edge]


class Subgraph(BaseModel):
    route_index: int
    route_node_labels: List[str]


class Routes(BaseModel):
    subgraphs: List[Subgraph]
    num_subgraphs: int


class SGPInput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    synth_graph_json: SynthGraphJsonInput = Field()
    unwanted_substances: List[str] = Field(default=[])
    unwanted_reactions: List[str] = Field(default=[])


class CustomConfig:
    @staticmethod
    def schema_extra(schema: Dict[str, Any], model: BaseModel) -> None:
        if "properties" in schema:
            schema["properties"]["synthesis_graph"] = {
                "title": "Synthesis Graph",
                "type": "object",
                "description": "A JSON representation of the synthesis graph",
            }


class TopNYieldSynthRoutesResult(ProfiledResult):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    search_params: TopNYieldSearch
    synth_graph: SynthGraph
    artifact_free_graph: SynthView
    combination_graphs: Optional[List[CombinationGraph]]
    route_candidates: Optional[SynthRouteContainer]
    routes: List[SynthRoute]
    time_info: TimeInfo
    summary: TopologyInfo
    availability: Optional[List[Availability]] = Field(default=None, description="Availability information for substances")

    def __init__(self, **data):
        super().__init__(**data)


class ReactionClassificationPrediction(BaseModel):
    rxsmiles: str
    predicted_classification: str
    probability: float


class ParseSynthGraphInput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    synthesis_graph_json: SynthGraphJson = Field()
    include_route_candidates: Optional[bool] = Field(default=False, description="Whether to return route candidates or not")
    include_combination_graphs: Optional[bool] = Field(default=False, description="Whether to return combination graphs or not")


class ParseSynthGraphOutput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    routes: List[SynthRoute] = Field(default=[], description="List of routes")
    route_candidates: Optional[List[SynthRoute]] = Field(default=None, description="List of route candidates")
    combination_graphs: Optional[List[CombinationGraph]] = Field(default=None, description="List of combination graphs")


class RouteLabelsInput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    aggregated_yield: Optional[float] = Field(default=None)
    predicted: bool = Field(default=False)
    source: Optional[str] = Field(default=None)
    route_index: Optional[int] = Field(default=None)
    route_status: Optional[str] = Field(default=None)
    method: Optional[str] = Field(default=None)
    route_node_labels: List[str] = Field(default=[], description="List of route node labels")


class AggregateYieldsInput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    repredict_yields: bool = Field(default=False, description="Whether to repredict yields or not")
    synthesis_graph_json: SynthGraphJson = Field(description="Synthesis graph input. Required to get yield information of reactions.")
    find_routes: bool = Field(default=False, description="Whether to find routes or not. If true, routes_json will be ignored.")
    routes: List[RouteLabelsInput] = Field(default=[], description="List of routes")


class AggregateYieldsOutput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    repredicted_yields: bool = Field(default=False, description="Whether yields were repredicted or not")
    synthesis_graph_json: SynthGraphJson = Field(description="Synthesis graph output.")
    found_routes: bool = Field(default=False, description="Whether routes were found or not")
    routes: List[SynthRoute] = Field(default=[], description="List of routes")