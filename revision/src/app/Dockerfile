FROM continuumio/miniconda3:latest

WORKDIR /workspace

# Install git so pip can install packages from GitHub URLs
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*

# Copy the conda environment file
COPY environment.yml /tmp/environment.yml

# Create the conda environment (syngps installed from GitHub)
RUN conda env create -f /tmp/environment.yml && conda clean --all -y

# Copy the app source
COPY . /workspace/app/

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "syngps-service", \
    "uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8002"]
