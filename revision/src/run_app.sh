#!/bin/bash

CONDA_ENV_FILE="./app/environment.yml"

# Check if the environment.yml file exists
if [ ! -f $CONDA_ENV_FILE ]; then
  echo "$CONDA_ENV_FILE file not found!"
  exit 1
fi

# Use awk to extract the name property from the environment.yml file
CONDA_ENV_NAME=$(awk '/^name:/ { print $2 }' $CONDA_ENV_FILE)

# Check if the name was found
if [ -z "$CONDA_ENV_NAME" ]; then
  echo "No name property found in $CONDA_ENV_FILE!"
  exit 1
else
  echo "Environment name: $CONDA_ENV_NAME"
fi

# Function to print usage information
usage() {
    echo "Usage: $0 [--dev] [--skip-env-setup] [--skip-wait-for-it] [-h|--help]"
    echo "  --dev                Run the FastAPI app in development mode with auto-reload"
    echo "  --skip-env-setup     Skip the Conda environment creation/update"
    echo "  --skip-wait-for-it   Skip waiting for the database to be available"
    echo "  -h, --help           Show this help message and exit"
}

# Function to load .env file
load_env_file() {
    if [ -f ./app/.env ]; then
        export $(grep -v '^#' ./app/.env | xargs)
    else
        echo "'.env' file not found. Using default environment variables."
    fi
}

# Function to convert DEBUG value to boolean
is_debug_mode() {
    # Convert DEBUG to lowercase
    debug_value=$(echo "$DEBUG" | tr '[:upper:]' '[:lower:]')
    if [ "$debug_value" = "true" ]; then
        echo true
    else
        echo false
    fi
}

# Function to check for conda or micromamba
check_conda_or_micromamba() {
    if command -v conda &> /dev/null; then
        echo "Conda is installed"
        CONDA_COMMAND="conda"
        EXTRA_RUN_OPTIONS="--no-capture-output"
    elif command -v micromamba &> /dev/null; then
        echo "Micromamba is installed"
        CONDA_COMMAND="micromamba"
        EXTRA_RUN_OPTIONS=""
    else
        echo "Neither Conda nor Micromamba is installed. Please install one of them."
        exit 1
    fi
    echo "CONDA_COMMAND: $CONDA_COMMAND"
}

# Function to create Conda environment if it does not exist
setup_conda_environment() {
    if $CONDA_COMMAND env list | grep -q "$CONDA_ENV_NAME"; then
        echo "Conda environment '$CONDA_ENV_NAME' already exists. Updating..."
        $CONDA_COMMAND env update -f "$CONDA_ENV_FILE"
    else
        echo "Creating Conda environment from $CONDA_ENV_FILE..."
        $CONDA_COMMAND env create -f "$CONDA_ENV_FILE"
    fi
}

# Function to start FastAPI app with Uvicorn using conda run
start_fastapi() {
    echo "Starting FastAPI app..."
    local mode=$1
    local debug_flag=$2

    if [ "$mode" = "development" ]; then
        exec $CONDA_COMMAND run $EXTRA_RUN_OPTIONS -n "$CONDA_ENV_NAME" uvicorn app.main:app --reload --host $APP_HOST --port $APP_PORT $debug_flag
    else
        exec $CONDA_COMMAND run $EXTRA_RUN_OPTIONS -n "$CONDA_ENV_NAME" uvicorn app.main:app --host $APP_HOST --port $APP_PORT $debug_flag
    fi
}

# Function to check if a port is open
wait_for_port() {
    local host=$1
    local port=$2
    local retries=${WAIT_FOR_IT_RETRIES:-30}
    local wait=${WAIT_FOR_IT_SLEEP:-2}

    echo "Polling $host:$port to be available, will retry $retries times with a $wait second wait..."
    for ((i=0; i<retries; i++)); do
        nc -z "$host" "$port" 2>/dev/null && return 0
        sleep "$wait"
    done

    echo "Timed out waiting for $host:$port after $((retries*wait)) seconds"
    return 1

}

# Parse arguments
mode="production"
skip_env_setup=false
skip_wait_for_it=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev) mode="development";;
        --skip-env-setup) skip_env_setup=true;;
        --skip-wait-for-it) skip_wait_for_it=true;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1;;
    esac
    shift
done

# Load the .env file
load_env_file

# Determine if debug mode is enabled
debug_mode=$(is_debug_mode)
if [ "$debug_mode" = true ]; then
    debug_flag="--log-level debug"
else
    debug_flag=""
fi

# Check if conda or micromamba is installed
check_conda_or_micromamba

# Setup Conda environment if not skipped
if [ "$skip_env_setup" = false ]; then
    echo "Setting up Conda environment..."
    setup_conda_environment
fi

if [ "$skip_wait_for_it" = false ]; then
    # Wait for MongoDB to be ready
    mongo_db_host="$MONGO_DB_HOST"
    mongo_db_port="$MONGO_DB_PORT"
    echo "Waiting for MongoDB to be ready..."
    wait_for_port "$mongo_db_host" "$mongo_db_port" || exit 1
    echo "MongoDB is ready..."

    # If GRAPHDB_BACKEND is not set default to 'memgraph'
    if [ -z "$GRAPHDB_BACKEND" ]; then
        export GRAPHDB_BACKEND="memgraph"
    fi

    # Check the graph database backend
    if [ "$GRAPHDB_BACKEND" = "neo4j" ]; then
        # Wait for Neo4j to be ready
        neo4j_host="$NEO4J_HOST"
        neo4j_port="$NEO4J_PORT"
        echo "Waiting for Neo4j to be ready..."
        wait_for_port "$neo4j_host" "$neo4j_port" || exit 1
        echo "Neo4j is ready..."
    elif [ "$GRAPHDB_BACKEND" = "memgraph" ]; then
        # Wait for Memgraph to be ready
        memgraph_host="$MEMGRAPH_HOST"
        memgraph_port="$MEMGRAPH_PORT"
        echo "Waiting for Memgraph to be ready..."
        wait_for_port "$memgraph_host" "$memgraph_port" || exit 1
        echo "Memgraph is ready..."
    else
        echo "Invalid value for GRAPHDB_BACKEND. Please set it to 'neo4j' or 'memgraph'."
        exit 1
    fi
else
    echo "Skipping waiting for MongoDB and graph database..."
fi
# Start the FastAPI app
start_fastapi $mode "$debug_flag"
