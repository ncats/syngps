import math
from pstats import Stats
from typing import Any, Dict, List, Optional, Union

from networkx import DiGraph
from pydantic import BaseModel, ConfigDict, Field, computed_field, model_validator
from src.aicplib.utils import (
    graph_to_json,
    graph_to_json_edges,
    graph_to_json_nodes,
    hash_graph,
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


class SynthGraphSearch(BaseModel):
    target_molecule_inchikey: str = Field(examples=["YTIQRXMAVJLXHH-UHFFFAOYSA-N"], min_length=1)
    reaction_steps: int = Field(default=2, examples=[2], description="Number of reaction steps to search away from target molecule")
    query_type: str = Field(
        default="shortest_path",
        examples=["shortest_path"],
        description='Type of query to be used for extracting graph data, options are "shortest_path". "shortest_path" (the default) will use shortest path to truncate much of the graph searching to enable faster query time.',
    )
    leaves_as_sm: bool = Field(default=True, examples=[True], description="Whether to consider leaf substance nodes as starting materials")
    include_graph_svgs: bool = Field(default=False, examples=[False], description="Whether to include SVGs of the graph in the response")
    graph_backend: str = Field(default="memgraph", examples=["memgraph", "neo4j"], description="Whether to use memgraph or neo4j backend")

    @property
    def search_depth(self) -> int:
        # Automatically calculate search_depth as 2 * reaction_steps
        return 2 * self.reaction_steps


class TopNYieldSearch(SynthGraphSearch):
    top_n_routes: int = Field(default=3, examples=[3])
    include_route_candidates: bool = Field(default=False, examples=[False], description="Whether to return route candidates or not")
    include_combination_graphs: bool = Field(default=False, examples=[False], description="Whether to return combination graphs or not")
    synthesis_graph_json: Optional[Dict[str, Any]] = Field(default=None, description="Optional synthesis graph input")


class SynthRoute(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    route_candidate: DiGraph = Field(exclude=True)
    aggregate_yield: Optional[float] = Field(default=None)
    predicted: bool = Field(default=False)
    source: Optional[str] = Field(default=None)
    route_index: int = Field(default=None)
    route_status: str = Field(default=None)
    method: str = Field(default=None)

    @computed_field
    def route_node_labels(self) -> List[str]:
        """Extract and return the node_label of all nodes in the route_candidate graph."""
        return [node_data.get("node_label") for _, node_data in self.route_candidate.nodes(data=True)]


class SynthRoutePreDefined(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    route_candidate: None = None
    aggregate_yield: Optional[float] = Field(default=None)
    predicted: bool = Field(default=False)
    source: Optional[str] = Field(default=None)
    route_index: int = Field(default=None)
    route_status: str = Field(default=None)
    method: str = Field(default=None)
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


class TopologyInfo(BaseModel):
    num_routes: int
    num_combination_graphs: int
    num_route_candidates: int


class SynthGraphJsonInput(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="ignore")

    nodes: List[Dict] = Field(default=None)
    edges: List[Dict] = Field(default=None)


class RouteAssemblyType(BaseModel):
    is_predicted: bool
    is_evidence: bool


class Node(BaseModel):
    node_label: str
    node_type: str
    uuid: str


class SubstanceNode(Node):
    node_type: str = "substance"
    inchikey: str
    canonical_smiles: str
    srole: str
    route_assembly_type: RouteAssemblyType


class Provenance(BaseModel):
    is_in_aicp: bool
    is_in_uspto_full: bool
    is_in_savi_130k: bool


class Validation(BaseModel):
    is_balanced: bool


class YieldInfo(BaseModel):
    yield_predicted: float
    yield_score: int


class ReactionNode(Node):
    node_type: str = "reaction"
    rxid: str
    rxsmiles: str
    yield_info: YieldInfo
    provenance: Provenance  # Contains is_in_aicp
    validation: Validation  # Contains is_balanced
    route_assembly_type: RouteAssemblyType


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


class AggregateYieldsInput(BaseModel):
    model_config = {"extra": "ignore"}
    repredict_yields: bool = Field(default=False)
    synth_graph: SynthGraphJson
    routes: Routes


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


class TimeInfo(BaseModel):
    synth_graph_time: Optional[float]
    enumerate_routes_time: Optional[float]
    aggregate_yield_time: Optional[float]
    rank_routes_time: Optional[float]


class TopNYieldSynthRoutesResult(ProfiledResult):
    # Define model configuration
    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    # Define model fields
    search_params: TopNYieldSearch
    synth_graph: SynthGraph
    artifact_free_graph: SynthView
    combination_graphs: SynthRouteContainer
    route_candidates: SynthRouteContainer
    routes: SynthRouteContainer
    time_info: TimeInfo
    summary: TopologyInfo

    def __init__(self, **data):
        super().__init__(**data)


class ReactionClassificationPrediction(BaseModel):
    rxsmiles: str
    predicted_classification: str
    probability: float
