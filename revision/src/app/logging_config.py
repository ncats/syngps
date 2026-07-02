import json
import logging
import os
from contextvars import ContextVar
from logging.config import dictConfig

# Define the transaction ID context variable
TID_CONTEXT_VAR = ContextVar("AICP-TID", default="N/A")


class TransactionIDFilter(logging.Filter):
    def filter(self, record):
        transaction_id = TID_CONTEXT_VAR.get()
        record.transaction_id = transaction_id if transaction_id is not None else "N/A"
        return True


def log_record_factory(*args, **kwargs):
    record = old_factory(*args, **kwargs)
    record.transaction_id = TID_CONTEXT_VAR.get("N/A")
    return record


# Store the original factory
old_factory = logging.getLogRecordFactory()
# Set the new factory
logging.setLogRecordFactory(log_record_factory)

# Load log level from environment variable
log_level = os.getenv("LOG_LEVEL", "INFO").upper()
debug_mode = os.getenv("DEBUG_MODE", "FALSE").lower() == "true"

# App logger name
LOGGING_NAME = "synthplanning_logger"

# Define the logging configuration dictionary
LOGGING_CONFIG = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "format": "[%(asctime)s] %(levelname)s: (TID: %(transaction_id)s) - %(message)s",
        },
        "uvicorn": {  # Adding a formatter for uvicorn logs
            "format": "[%(asctime)s] %(levelname)s: %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "default",
        },
    },
    "loggers": {
        LOGGING_NAME: {
            "handlers": ["console"],
            "level": log_level,
            "propagate": False,
        },
        "uvicorn": {  # Applying the formatter to uvicorn logs
            "handlers": ["console"],
            "level": log_level,
            "propagate": False,
            "formatter": "uvicorn",
        },
        "uvicorn.error": {
            "handlers": ["console"],
            "level": log_level,
            "propagate": False,
            "formatter": "uvicorn",
        },
        "uvicorn.access": {
            "handlers": ["console"],
            "level": log_level,
            "propagate": False,
            "formatter": "uvicorn",
        },
        "aicplib": {
            "level": os.getenv("AICPLIB_LOG_LEVEL", log_level),
            "handlers": ["console"],
            "propagate": True,
            "qualname": "aicplib",
        }
    },
}


def setup_logging():
    # Apply the logging configuration
    dictConfig(LOGGING_CONFIG)

    # Create a logger for the current module
    logger = logging.getLogger(__name__)
    logger.addFilter(TransactionIDFilter())
    logger.setLevel(log_level)  # Ensure the logger level is set to the appropriate level

    # Add TransactionIDFilter to uvicorn loggers
    for logger_name in ["uvicorn", "uvicorn.error", "uvicorn.access"]:
        uvicorn_logger = logging.getLogger(logger_name)
        uvicorn_logger.addFilter(TransactionIDFilter())

    # Add filter & level to aicplib logger
    aicplib_logger = logging.getLogger("aicplib")
    aicplib_logger.setLevel(os.getenv("AICPLIB_LOG_LEVEL", log_level))
    aicplib_logger.addFilter(TransactionIDFilter())
    aicplib_logger.debug("aicplib logger initialized (level=%s)", aicplib_logger.level)

# Setup logging
setup_logging()

# Test the logging
logger = logging.getLogger(LOGGING_NAME)

# Log the logging configuration if in debug mode\
logger.debug(f"Logging configuration: \n{json.dumps(LOGGING_CONFIG, indent=2)}")
