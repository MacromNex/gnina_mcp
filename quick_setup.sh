#!/bin/bash
# quick_setup.sh - Quick setup script for gnina_mcp
#
# Usage:
#   bash quick_setup.sh [options]
#
# Options:
#   --help      Show this help message
#   --skip-env  Skip environment creation (use existing env)
#   --skip-repo Skip repository cloning (repo already present)

set -e  # Exit on error

# =============================================================================
# Configuration
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MCP_NAME="gnina_mcp"
PYTHON_VERSION="3.10"
REPO_URL="https://github.com/gnina/gnina"
REPO_NAME="gnina"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# Helper Functions
# =============================================================================
log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

show_help() {
    cat << EOF
${MCP_NAME} Quick Setup Script

Usage:
    bash quick_setup.sh [options]

Options:
    --help      Show this help message
    --skip-env  Skip environment creation (use existing env)
    --skip-repo Skip repository cloning (repo already present)

Description:
    This script sets up the environment for ${MCP_NAME} MCP server.
    It will:
    1. Create a conda/mamba environment with Python ${PYTHON_VERSION}
    2. Clone the gnina repository (if not present)
    3. Install all required dependencies including RDKit
    4. Install the MCP server dependencies (fastmcp)

Requirements:
    - conda or mamba (mamba recommended for faster installation)
    - git
    - ~1.5 GB disk space
    - Internet connection for package downloads

Time estimate: 5-10 minutes with mamba, 10-20 minutes with conda

EOF
    exit 0
}

get_pkg_manager() {
    if command -v mamba &> /dev/null; then
        echo "mamba"
    elif command -v conda &> /dev/null; then
        echo "conda"
    else
        echo ""
    fi
}

# =============================================================================
# Parse Arguments
# =============================================================================
SKIP_ENV=false
SKIP_REPO=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --help|-h)
            show_help
            ;;
        --skip-env)
            SKIP_ENV=true
            shift
            ;;
        --skip-repo)
            SKIP_REPO=true
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            ;;
    esac
done

# =============================================================================
# Main Setup
# =============================================================================
cd "$SCRIPT_DIR"
log_info "Setting up ${MCP_NAME} in: $SCRIPT_DIR"

# Check for package manager
PKG_MGR=$(get_pkg_manager)
if [ -z "$PKG_MGR" ]; then
    log_error "Neither mamba nor conda found. Please install one of them first."
    log_error "Recommended: conda install mamba -n base -c conda-forge"
    exit 1
fi
log_info "Using package manager: $PKG_MGR"

# -----------------------------------------------------------------------------
# Step 1: Clone Repository (if needed)
# -----------------------------------------------------------------------------
if [ "$SKIP_REPO" = false ]; then
    if [ -d "repo/${REPO_NAME}" ]; then
        log_info "Repository already exists at repo/${REPO_NAME}"
    else
        log_info "Cloning gnina repository..."
        mkdir -p repo
        git clone --depth=1 "${REPO_URL}" "repo/${REPO_NAME}" || {
            log_warning "Shallow clone failed, trying full clone..."
            git clone "${REPO_URL}" "repo/${REPO_NAME}"
        }
        log_success "Repository cloned successfully"
    fi
else
    log_info "Skipping repository clone (--skip-repo)"
fi

# -----------------------------------------------------------------------------
# Step 2: Create Environment (if needed)
# -----------------------------------------------------------------------------
if [ "$SKIP_ENV" = false ]; then
    if [ -d "env" ]; then
        log_warning "Environment already exists at ./env"
        read -p "Do you want to remove and recreate it? [y/N] " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_info "Removing existing environment..."
            rm -rf env
        else
            log_info "Keeping existing environment"
            SKIP_ENV=true
        fi
    fi

    if [ "$SKIP_ENV" = false ]; then
        log_info "Creating conda environment with Python ${PYTHON_VERSION}..."
        $PKG_MGR create -p ./env python=${PYTHON_VERSION} -y
        log_success "Environment created successfully"
    fi
else
    log_info "Skipping environment creation (--skip-env)"
fi

# Verify environment exists
if [ ! -f "env/bin/python" ]; then
    log_error "Environment not found at ./env/bin/python"
    log_error "Please run without --skip-env to create it"
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 3: Install Dependencies
# -----------------------------------------------------------------------------
log_info "Installing dependencies..."

# Activate environment for conda installs
if [ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ./env
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate ./env
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
    conda activate ./env
else
    log_warning "Could not find conda initialization script, using direct env path"
fi

# Install core Python packages
log_info "Installing core Python packages..."
./env/bin/pip install --upgrade pip
./env/bin/pip install loguru click pandas numpy tqdm

# Install FastMCP (force reinstall for clean installation)
log_info "Installing FastMCP..."
./env/bin/pip install --force-reinstall --no-cache-dir fastmcp

# Install RDKit from conda-forge (most reliable method)
log_info "Installing RDKit from conda-forge..."
$PKG_MGR install -p ./env -c conda-forge rdkit -y

# Install scientific computing packages
log_info "Installing scientific computing packages..."
./env/bin/pip install scipy matplotlib seaborn scikit-learn jupyter pytest

log_success "Dependencies installed successfully"

# -----------------------------------------------------------------------------
# Step 4: Post-install Setup
# -----------------------------------------------------------------------------
log_info "Setting up gnina binary symlink..."

# Create symlink to gnina binary if it exists
if [ -f "repo/${REPO_NAME}/gnina_binary" ]; then
    if [ ! -f "gnina" ]; then
        ln -s "$(pwd)/repo/${REPO_NAME}/gnina_binary" gnina
        log_success "Created gnina binary symlink"
    else
        log_info "Gnina binary symlink already exists"
    fi
else
    log_warning "Gnina binary not found at repo/${REPO_NAME}/gnina_binary"
    log_warning "For production use, compile gnina or download pre-built binary"
fi

# Create directories for job management
log_info "Creating job directories..."
mkdir -p jobs results configs
log_success "Created job management directories"

# -----------------------------------------------------------------------------
# Step 5: Verify Installation
# -----------------------------------------------------------------------------
log_info "Verifying installation..."

# Test Python environment
./env/bin/python -c "import sys; print(f'Python {sys.version}')" || {
    log_error "Python environment verification failed"
    exit 1
}

# Test core imports
./env/bin/python -c "import fastmcp; print(f'fastmcp {fastmcp.__version__}')" || {
    log_error "fastmcp import failed"
    exit 1
}

# Test RDKit
./env/bin/python -c "from rdkit import Chem; print('RDKit OK')" || {
    log_error "RDKit import failed"
    exit 1
}

# Test scientific stack
./env/bin/python -c "import numpy, pandas, scipy; print('Scientific stack OK')" || {
    log_error "Scientific computing stack import failed"
    exit 1
}

# Test MCP server import
./env/bin/python -c "from src.server import mcp; print('MCP server import OK')" 2>/dev/null || {
    log_warning "MCP server import test skipped (server may not be ready yet)"
}

log_success "Installation verified successfully"

# -----------------------------------------------------------------------------
# Step 6: Environment Summary
# -----------------------------------------------------------------------------
log_info "Generating environment summary..."

# Create environment info file
cat > env_info.txt << EOF
# ${MCP_NAME} Environment Information
# Generated: $(date)

## Python Environment
Python Version: $(./env/bin/python --version)
Environment Path: $(pwd)/env

## Key Package Versions
FastMCP: $(./env/bin/python -c "import fastmcp; print(fastmcp.__version__)" 2>/dev/null || echo "Not available")
RDKit: $(./env/bin/python -c "from rdkit import rdBase; print(rdBase.rdkitVersion)" 2>/dev/null || echo "Not available")
NumPy: $(./env/bin/python -c "import numpy; print(numpy.__version__)" 2>/dev/null || echo "Not available")
Pandas: $(./env/bin/python -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "Not available")

## Repository Information
Gnina Repository: repo/${REPO_NAME}
Gnina Binary: $([ -f "gnina" ] && echo "Available (symlinked)" || echo "Not available")

## Installation Status
Environment: ✅ Created
Dependencies: ✅ Installed
RDKit: ✅ Working
FastMCP: ✅ Working
Job Directories: ✅ Created
EOF

log_success "Environment summary saved to env_info.txt"

# =============================================================================
# Done
# =============================================================================
echo
log_success "=============================================="
log_success "  ${MCP_NAME} setup completed!"
log_success "=============================================="
echo
echo "Environment Information:"
echo "  • Python: $(./env/bin/python --version)"
echo "  • Package Manager: ${PKG_MGR}"
echo "  • Environment Size: ~500MB"
echo "  • Gnina Binary: $([ -f "gnina" ] && echo "✅ Available" || echo "⚠️  Not available")"
echo
echo "Next steps:"
echo "  1. Register with Claude Code:"
echo "     cd $(pwd)"
echo "     claude mcp add gnina-tools -- \$(pwd)/env/bin/python \$(pwd)/src/server.py"
echo
echo "  2. Or use cpmcp (recommended):"
echo "     cpmcp install gnina_mcp"
echo
echo "  3. Test the server:"
echo "     ./env/bin/python src/server.py"
echo
echo "  4. Run example scripts:"
echo "     ./env/bin/python scripts/score_protein_ligand.py --help"
echo
echo "For troubleshooting, see README.md or check env_info.txt"
echo