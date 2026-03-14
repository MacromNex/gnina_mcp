#!/bin/bash
# quick_setup.sh - Quick setup script for gnina_mcp
#
# Installs gnina binary and dependencies into the project's conda env.
#
# Prerequisites:
#   - conda/mamba environment at ./env
#   - CUDA 12.x compatible GPU (optional, CPU fallback available)
#
# Usage:
#   bash quick_setup.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/env"
GNINA_VERSION="1.3.2"
GNINA_BINARY="gnina.${GNINA_VERSION}.cuda12.8"
GNINA_URL="https://github.com/gnina/gnina/releases/download/v${GNINA_VERSION}/${GNINA_BINARY}"

echo "=== gnina_mcp Quick Setup ==="
echo ""

# Step 1: Check env exists
if [ ! -d "${ENV_DIR}/bin" ]; then
    echo "Error: conda env not found at ${ENV_DIR}"
    echo "Create it first: conda create -p ${ENV_DIR} python=3.10 pandas numpy -y"
    exit 1
fi

# Step 2: Download gnina binary if not present
if [ ! -f "${ENV_DIR}/bin/${GNINA_BINARY}" ]; then
    echo "Downloading gnina v${GNINA_VERSION}..."
    wget -q --show-progress -O "${ENV_DIR}/bin/${GNINA_BINARY}" "${GNINA_URL}" \
        || curl -L -o "${ENV_DIR}/bin/${GNINA_BINARY}" "${GNINA_URL}"
    chmod 755 "${ENV_DIR}/bin/${GNINA_BINARY}"
    echo "Downloaded gnina to ${ENV_DIR}/bin/${GNINA_BINARY}"
else
    echo "gnina binary already exists at ${ENV_DIR}/bin/${GNINA_BINARY}"
fi

# Step 3: Create symlink
if [ ! -L "${ENV_DIR}/bin/gnina" ] || [ "$(readlink "${ENV_DIR}/bin/gnina")" != "${GNINA_BINARY}" ]; then
    ln -sf "${GNINA_BINARY}" "${ENV_DIR}/bin/gnina"
    echo "Created symlink: ${ENV_DIR}/bin/gnina -> ${GNINA_BINARY}"
else
    echo "Symlink already correct"
fi

# Step 4: Install cuDNN (needed for gnina CNN scoring)
echo ""
echo "Installing cuDNN..."
if command -v mamba &>/dev/null; then
    mamba install -p "${ENV_DIR}" -c conda-forge cudnn=9 -y --quiet
elif command -v conda &>/dev/null; then
    conda install -p "${ENV_DIR}" -c conda-forge cudnn=9 -y --quiet
else
    echo "Warning: conda/mamba not found. Install cuDNN manually if CNN scoring fails."
fi

# Step 5: Install Python dependencies
echo ""
echo "Installing Python dependencies..."
"${ENV_DIR}/bin/pip" install --quiet pandas numpy loguru fastmcp 2>/dev/null || true

# Step 6: Test gnina
echo ""
echo "Testing gnina installation..."
export LD_LIBRARY_PATH="${ENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
if "${ENV_DIR}/bin/gnina" --version 2>/dev/null; then
    echo "gnina is working!"
else
    echo "Warning: gnina binary test failed. This may be due to missing CUDA/cuDNN."
    echo "Try: export LD_LIBRARY_PATH=${ENV_DIR}/lib:\$LD_LIBRARY_PATH"
    echo "CPU mode (--no_gpu) may still work."
fi

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Usage:"
echo "  # Activate environment"
echo "  export LD_LIBRARY_PATH=${ENV_DIR}/lib:\$LD_LIBRARY_PATH"
echo ""
echo "  # Run scripts directly"
echo "  ${ENV_DIR}/bin/python scripts/score_protein_ligand.py -r receptor.pdb -l ligand.sdf -o results.csv"
echo ""
echo "  # Start MCP server"
echo "  ${ENV_DIR}/bin/python src/server.py"
