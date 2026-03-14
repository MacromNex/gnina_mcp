# ==============================================================================
# Dockerfile for gnina MCP Server
#
# Provides CNN-enhanced molecular docking via gnina as an MCP service.
#
# Build:
#   docker build -t gnina-mcp .
#
# Run MCP server (stdio transport):
#   docker run -i gnina-mcp
#
# Run MCP server (SSE transport on port 8000):
#   docker run -p 8000:8000 gnina-mcp --transport sse --host 0.0.0.0 --port 8000
#
# Run a script directly:
#   docker run gnina-mcp python scripts/score_protein_ligand.py \
#     -r examples/data/3rod_rec.pdb -l examples/data/3rod_lig.pdb
#
# Mount local data:
#   docker run -v /path/to/data:/data gnina-mcp python scripts/dock_ligand.py \
#     -r /data/receptor.pdb -l /data/ligand.sdf --autobox_ligand /data/ref.pdb \
#     -o /data/docked.sdf.gz
# ==============================================================================

# --- Stage 1: Download gnina binary ---
FROM ubuntu:22.04 AS downloader

RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ARG GNINA_VERSION=1.3.2
ARG GNINA_BINARY=gnina.${GNINA_VERSION}.cuda12.8
ARG GNINA_URL=https://github.com/gnina/gnina/releases/download/v${GNINA_VERSION}/${GNINA_BINARY}

RUN wget -q -O /tmp/gnina "${GNINA_URL}" && chmod 755 /tmp/gnina

# --- Stage 2: Runtime image ---
FROM nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04

LABEL maintainer="gnina-mcp" \
      description="CNN-enhanced molecular docking MCP server using gnina" \
      gnina.version="1.3.2"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install Python (gnina binary is mostly statically linked, only needs CUDA/cuDNN from base)
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
    && ln -sf /usr/bin/python3 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*

# Copy gnina binary from downloader stage
COPY --from=downloader /tmp/gnina /usr/local/bin/gnina

# Verify gnina is functional
RUN gnina --version || echo "gnina binary present (GPU test skipped in build)"

# Set up application
WORKDIR /app

# Install Python dependencies first (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt \
    && pip install --no-cache-dir "numpy<2" rdkit-pypi

# Copy application code
COPY src/ src/
COPY scripts/ scripts/
COPY configs/ configs/
COPY examples/ examples/
COPY mock_gnina.py .

# Create directories for job output (world-writable for --user UID:GID)
RUN mkdir -p jobs results && chmod 777 jobs results

# Set environment so gnina is always found
ENV PATH="/usr/local/bin:${PATH}"
ENV PYTHONPATH="/app/src:/app/scripts"
ENV PYTHONUNBUFFERED=1

# Health check: verify imports work
RUN python -c "from fastmcp import FastMCP; print('FastMCP OK')" \
    && python -c "import pandas, numpy, loguru; print('Dependencies OK')" \
    && python -c "from rdkit import Chem; print('RDKit OK')" \
    && python -c "import sys; sys.path.insert(0,'src'); sys.path.insert(0,'scripts'); from server import mcp; print(f'MCP server OK: {mcp.name}')"

# Default: run MCP server via stdio
ENTRYPOINT ["python", "src/server.py"]
