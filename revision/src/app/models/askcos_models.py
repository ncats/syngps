import re
import uuid
from typing import Any, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

################################################
# CUSTOM MODELS TO HANDLE ASKCOS API RESPONSES #
################################################


class ModuleStatus(BaseModel):
    name: str
    description: str
    available_model_names: List[str]
    ready: bool


class ModuleStatusesModel(BaseModel):
    modules: List[ModuleStatus]
    all_healthy: bool


class AskcosAtommapResult(BaseModel):
    mapped_rxn: str
    confidence: float


class AskcosAtommapResponse(BaseModel):
    status_code: int
    message: str
    result: List[AskcosAtommapResult]


class AskcosTreeRequest(BaseModel):
    model_config = ConfigDict(extra="ignore")
    target_molecule_smiles: str = Field(description="Target molecule SMILES for tree search", examples=["CN(C)CCOC(c1ccccc1)c1ccccc1"])
    max_tree_depth: int = Field(default=5, description="Maximum tree depth for tree search", examples=[5])
    max_trees: int = Field(default=25, description="Maximum number of trees to search", examples=[25])


################################
# MODELS FROM ASKCOS CODE BASE #
################################
# Source: https://github.com/ncats/aspire_askcos_askcos2_core/tree/main
# Note: In ASKCOS code the following models were split between multiple files, joined together here.


def to_pascal(snake: str) -> str:
    """Convert a snake_case string to PascalCase.

    Args:
        snake: The string to convert.

    Returns:
        The PascalCase string.
    """
    camel = snake.title()
    return re.sub('([0-9A-Za-z])_(?=[0-9A-Z])', lambda m: m.group(1), camel)


def to_camel(snake: str) -> str:
    """Convert a snake_case string to camelCase.

    Args:
        snake: The string to convert.

    Returns:
        The converted camelCase string.
    """
    camel = to_pascal(snake)
    return re.sub('(^_*[A-Z])', lambda m: m.group(1).lower(), camel)


class LowerCamelAliasModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def populate_by_alias_if_needed(cls, data: dict[str, Any]) -> dict[str, Any]:
        for k, field_info in cls.model_fields.items():
            if k not in data and to_camel(k) in data:
                data[k] = data.pop(to_camel(k))
        return data


class BaseResponse(BaseModel):
    status_code: int
    message: str
    result: Any


class RetroBackendOption(BaseModel):
    retro_backend: Literal["augmented_transformer", "graph2smiles", "template_relevance", "retrosim"] = Field(
        default="template_relevance", description="backend for one-step retrosynthesis"
    )
    retro_model_name: str = Field(default="reaxys", description="backend model name for one-step retrosynthesis")
    max_num_templates: int = Field(default=50, description="number of templates to consider")
    max_cum_prob: float = Field(default=0.995, description="maximum cumulative probability of templates")
    attribute_filter: list[dict[str, Any]] = Field(
        default_factory=list, description="template attribute filter to apply before template application", examples=[[]]
    )
    threshold: float = Field(default=0.3, description="threshold for similarity; " "used for retrosim only")
    top_k: int = Field(default=10, description="filter for k results returned; " "used for retrosim only")

    @field_validator("retro_model_name", mode="after")
    @classmethod
    def validate_retro_model_name(cls, value, info):
        backend = info.data.get("retro_backend")
        if backend is None:
            raise ValueError("retro_backend not supplied!")

        valid_model_names = {
            "template_relevance": ["bkms_metabolic", "cas", "pistachio", "pistachio_ringbreaker", "reaxys", "reaxys_biocatalysis"],
            "augmented_transformer": ["cas", "pistachio_23Q3", "USPTO_FULL"],
            "graph2smiles": ["cas", "pistachio_23Q3", "USPTO_FULL"],
            "retrosim": ["USPTO_FULL", "bkms"],
        }

        if backend in valid_model_names and value not in valid_model_names[backend]:
            raise ValueError(f"Unsupported retro_model_name {value} for {backend}")
        return value


class ClusterSetting(LowerCamelAliasModel):
    feature: Literal["original", "outcomes", "all"] = "original"
    cluster_method: Literal["rxn_class", "hdbscan", "kmeans"] = "hdbscan"
    fp_type: Literal["morgan"] = Field(default="morgan", description="fingerprint type")
    fp_length: int = Field(default=512, description="fingerprint bits")
    fp_radius: int = Field(default=1, description="fingerprint radius")
    classification_threshold: float = Field(default=0.2, description="threshold to classify a reaction as unknown when " "clustering by 'rxn_class'")


class ExpandOneOptions(LowerCamelAliasModel):
    # aliasing to v1 fields
    template_max_count: int = Field(default=100, alias="template_count")
    template_max_cum_prob: float = Field(default=0.995, alias="max_cum_template_prob")
    banned_chemicals: list[str] = Field(
        default_factory=list, description="banned chemicals (in addition to user banned chemicals)", alias="forbidden_molecules", examples=[[]]
    )
    banned_reactions: list[str] = Field(
        default_factory=list, description="banned reactions (in addition to user banned reactions)", alias="known_bad_reactions", examples=[[]]
    )

    retro_backend_options: list[RetroBackendOption] = Field(default=[RetroBackendOption()], description="list of retro strategies to run in series")
    use_fast_filter: bool = Field(default=True, description="whether to filter the results with the fast filter")
    filter_threshold: float = Field(default=0.75, description="threshold for the fast filter")
    retro_rerank_backend: Literal["relevance_heuristic", "scscore"] | None = Field(default=None, description="backend for retro rerank")
    cluster_precursors: bool = Field(default=False, description="whether to cluster proposed precursors")
    cluster_setting: ClusterSetting = Field(default_factory=ClusterSetting, description="settings for clustering")
    extract_template: bool = Field(default=False, description="whether to extract templates " "(mostly for template-free suggestions)")
    return_reacting_atoms: bool = Field(default=True, description="whether to return the indices of reacting atoms")
    selectivity_check: bool = Field(default=False, description="whether to perform quick selectivity check " "by reverse application of the forward template")

    class Config:
        populate_by_name = True


class BuildTreeOptions(LowerCamelAliasModel):
    expansion_time: int = Field(default=30, description="max time for tree search in seconds", examples=[10])
    max_iterations: int | None = Field(default=None, description="max number of iterations", examples=[None])
    max_chemicals: int | None = Field(default=None, description="max number of chemicals to explore", examples=[None])
    max_reactions: int | None = Field(default=None, description="max number of reactions to explore", examples=[None])
    max_templates: int | None = Field(default=None, description="max number of templates to explore", examples=[None])
    max_branching: int = Field(default=25, description="max number of branching")
    max_depth: int = Field(default=5, description="max tree depth")
    exploration_weight: float = Field(default=1.0, description="weight for exploration (vs. exploitation)")
    return_first: bool = Field(default=False, description="whether to stop when the first buyable path is found")
    max_trees: int = Field(default=25, description="max number of buyable paths to explore")

    # a bunch of termination logic. These were passed directly from the front end.
    # Grouping happens at the Django side. Let's tentatively keep that pattern for now.
    buyable_logic: Literal["none", "and", "or"] | None = Field(default="and", description="logic type for buyable termination")
    max_ppg_logic: Literal["none", "and", "or"] | None = Field(default="none", description="logic type for price based termination")
    max_ppg: float | None = Field(default=None, description="maximum price for price based termination")
    max_scscore_logic: Literal["none", "and", "or"] | None = Field(default="none", description="logic type for synthetic complexity termination")
    max_scscore: float | None = Field(default=None, description="maximum scscore for synthetic complexity termination")
    chemical_property_logic: Literal["none", "and", "or"] | None = Field(default="none", description="logic type for chemical property termination")
    max_chemprop_c: int | None = Field(default=None, description="maximum carbon count for termination")
    max_chemprop_n: int | None = Field(default=None, description="maximum nitrogen count for termination")
    max_chemprop_o: int | None = Field(default=None, description="maximum oxygen count for termination")
    max_chemprop_h: int | None = Field(default=None, description="maximum hydrogen count for termination")
    chemical_popularity_logic: Literal["none", "and", "or"] | None = Field(default="none", description="logic type for chemical popularity termination")
    min_chempop_reactants: int | None = Field(default=5, description="minimum reactant precedents for termination")
    min_chempop_products: int | None = Field(default=5, description="minimum product precedents for termination")

    buyables_source: list[str] | None = Field(default=None, description="list of source(s) to consider when looking up buyables", examples=[None])
    custom_buyables: list[str] | None = Field(default=None, description="list of chemicals to consider as buyable", examples=[[]])

    use_value_network: bool | None = Field(default=False, description="whether to use value_network for Vm computation")


class EnumeratePathsOptions(LowerCamelAliasModel):
    path_format: Literal["json", "graph"] = "json"
    json_format: Literal["treedata", "nodelink"] = "nodelink"
    sorting_metric: Literal["plausibility", "number_of_starting_materials", "number_of_reactions", "score"] = "plausibility"
    validate_paths: bool = True
    score_trees: bool = False
    cluster_trees: bool = False
    cluster_method: Literal["hdbscan", "kmeans"] = "hdbscan"
    min_samples: int = 5
    min_cluster_size: int = 5
    paths_only: bool = False
    max_paths: int = 10


class TreeSearchInput(LowerCamelAliasModel):
    backend: Literal["mcts", "retro_star"] = Field(default="mcts", description="backend for tree search")
    smiles: str = Field(description="target SMILES for tree search", examples=["CN(C)CCOC(c1ccccc1)c1ccccc1"])
    description: str | None = Field(
        default="",
        description="description of the tree search task",
    )
    tags: str | None = Field(
        default="",
        description="tags of the tree search task",
    )
    expand_one_options: ExpandOneOptions = Field(default_factory=ExpandOneOptions, description="options for one-step expansion")
    build_tree_options: BuildTreeOptions = Field(default_factory=BuildTreeOptions, description="options for tree search")
    enumerate_paths_options: EnumeratePathsOptions = Field(
        default_factory=EnumeratePathsOptions, description="options for path enumeration once the tree is built"
    )
    run_async: bool = False
    result_id: str = str(uuid.uuid4())


class TreeSearchResult(BaseModel):
    stats: dict[str, Any] | None
    paths: Optional[list[dict[str, Any]]] | None = None
    graph: Optional[dict[str, Any]] | None = None
    version: int | str | None = 2
    result_id: str = ""


class TreeSearchOutput(BaseModel):
    placeholder: str


class TreeSearchResponse(BaseResponse):
    result: TreeSearchResult | None

class ContextPredictionRequest(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="allow", arbitrary_types_allowed=True)

    reactants: str
    products: str
    with_smiles: bool = False
    return_scores: bool = False
    num_results: int = 5


class ContextPredictionOption(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="allow", arbitrary_types_allowed=True)

    temperature: Optional[float]
    solvent: Optional[str]
    reagent: Optional[str]
    catalyst: Optional[str]
    solvent_score: Optional[float]
    best: Optional[bool]
    solvent_name_only: Optional[str]
    reagent_name_only: Optional[str]
    catalyst_name_only: Optional[str]
    score: Optional[float]

class ContextPredictionResponse(BaseModel):
    # Define model configuration
    model_config = ConfigDict(extra="allow", arbitrary_types_allowed=True)

    task_id: str
    state: str
    complete: bool
    failed: bool
    message: str | None = None
    output: List[ContextPredictionOption]

class AskcosConditionPredictRequest(BaseModel):
    rxsmiles: str = Field(default="", title="RXSMILES", description="RXSmiles to predict", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])
    num_results: int = 5

class AskcosConditionPredictResponse(BaseModel):
    rxsmiles: str = Field(default="", title="RXSMILES", description="RXSmiles to predict", examples=["CCO.CC(=O)O>>CC(=O)OCC.O"])
    condition_options: List[ContextPredictionOption]