import os
from pathlib import Path
from pydantic import Field, ValidationError, field_validator
from pydantic_settings import BaseSettings
from app.logging_config import logger


class AppConfig(BaseSettings):
    # Application settings
    app_port: int = Field(default=8002, alias="APP_PORT")
    app_host: str = Field(default="0.0.0.0", alias="APP_HOST")
    app_prefix: str = Field(default="", alias="APP_PREFIX")
    app_version: str = Field(alias="APP_VERSION")
    log_level: str = Field(default="INFO", alias="LOG_LEVEL")
    debug_mode: bool = Field(default=False, alias="DEBUG")
    smiles_encryption_key: str = Field(alias="SMILES_ENCRYPTION_KEY")
    models_dir: str = Field(alias="MODELS_DIR")

    # MongoDB settings
    mongo_db_host: str = Field(alias="MONGO_DB_HOST")
    mongo_db_port: int = Field(default=27017, alias="MONGO_DB_PORT")
    mongo_db_name: str = Field(alias="MONGO_DB_NAME")
    mongo_db_user: str | None = Field(default=None, alias="MONGO_DB_USER")
    mongo_db_pass: str | None = Field(default=None, alias="MONGO_DB_PASS")

    # Graph DB Backend
    graph_db_backend: str = Field(default="memgraph", alias="GRAPHDB_BACKEND")

    # Memgraph settings
    memgraph_host: str = Field(alias="MEMGRAPH_HOST")
    memgraph_port: int = Field(default=7687, alias="MEMGRAPH_PORT")
    memgraph_user: str | None = Field(default=None, alias="MEMGRAPH_USER")
    memgraph_pass: str | None = Field(default=None, alias="MEMGRAPH_PASS")
    memgraph_conn_encrypted: bool = Field(default=False, alias="MEMGRAPH_CONN_ENCRYPTED")

    # ASKCOS Connection Settings
    askcos_base_url: str = Field(default="https://askcos.mit.edu/api", alias="ASKCOS_BASE_URL")

    # Inventory file paths for adapters
    askcos_inv_file_path: str | None = Field(default=None, alias="ASKCOS_INV_FILE_PATH")
    custom_stock_inv_file_path: str | None = Field(default=None, alias="CUSTOM_STOCK_INV_FILE_PATH")

    # Validators can be added if needed
    @field_validator("log_level")
    def validate_log_level(cls, value):
        valid_levels = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
        if value.upper() not in valid_levels:
            raise ValueError(f"Invalid log level: {value}")
        return value.upper()

    # Validate app api prefix
    @field_validator("app_prefix")
    def validate_app_prefix(cls, value):
        if len(value) > 0 and not value.startswith("/"):
            raise ValueError(f"Invalid app prefix, must start with '/': {value}")
        return value

    @field_validator("askcos_base_url")
    def validate_askcos_base_url(cls, value: str) -> str:
        # Keep base URL stable for string concatenation with endpoints like "/admin/...".
        return value.strip().rstrip("/")

    @field_validator("models_dir")
    def validate_models_dir(cls, value):
        if not os.path.exists(value):
            raise ValueError(f"Models directory does not exist: {value}")
        if not os.path.isdir(value):
            raise ValueError(f"Models path is not a directory: {value}")
        if not os.access(value, os.R_OK):
            raise ValueError(f"Models directory is not readable: {value}")
        return os.path.abspath(value)

    @field_validator("debug_mode", "memgraph_conn_encrypted")
    def validate_boolean(cls, value):
        if isinstance(value, str):
            return value.lower() in ("true", "1", "t", "yes", "y")
        return value

    @property
    def mongo_connection_uri(self) -> str:
        if self.mongo_db_user and self.mongo_db_pass:
            return f"mongodb://{self.mongo_db_user}:{self.mongo_db_pass}@{self.mongo_db_host}:{self.mongo_db_port}"
        return f"mongodb://{self.mongo_db_host}:{self.mongo_db_port}"

    @property
    def memgraph_connection_uri(self) -> str:
        return f"bolt://{self.memgraph_host}:{self.memgraph_port}"
    
    class Config:
        extra = "allow"


# Determine if .env file exists
env_path = Path(".env")

try:
    if env_path.exists():
        class EnvAppConfig(AppConfig):
            class Config(AppConfig.Config):
                env_file = ".env"

        APP_CONFIG = EnvAppConfig() # type: ignore
        logger.info("Loaded configuration from .env file.")
    else:
        APP_CONFIG = AppConfig() # type: ignore
        logger.info("Loaded configuration from environment variables.")

    logger.debug(f"Configuration loaded: {APP_CONFIG.model_dump()}")
except ValidationError as e:
    logger.error(f"Configuration error: {e}")
    raise