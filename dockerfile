# ============================================================
# Immuno-Cancer Pipeline
# Author: Oguzhan Begik, PhD
# Description: Reproducible Snakemake-based TCGA Immune Survival Analysis
# ============================================================

FROM continuumio/miniconda3:latest

LABEL maintainer="Oguzhan Begik <oguzhanbegik@gmail.com>"
LABEL description="TCGA immune survival pipeline using Snakemake + Python"

# Set working directory
WORKDIR /app

# Copy environment definition
COPY environment.yml /app/environment.yml

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libxml2-dev libxslt-dev zlib1g-dev \
    libjpeg-dev libpng-dev libgl1 \
    poppler-utils \
    && rm -rf /var/lib/apt/lists/*

# Create the conda environment and install SerpAPI from GitHub (not PyPI)
RUN conda env create -f environment.yml && \
    conda run -n snakemake_env pip install --no-cache-dir torch transformers accelerate && \
    conda clean -afy

SHELL ["conda", "run", "-n", "snakemake_env", "/bin/bash", "-c"]

# Ensure snakemake is available in PATH
ENV PATH /opt/conda/envs/snakemake_env/bin:$PATH

# Copy all pipeline components
COPY Snakefile /app/
COPY scripts /app/scripts
COPY config /app/config
COPY assets /app/assets
COPY data /app/data
COPY results /app/results

# Default command: run snakemake
ENTRYPOINT ["conda", "run", "-n", "snakemake_env", "snakemake"]
CMD ["--cores", "2", "--snakefile", "Snakefile", "--configfile", "config/config.yaml"]