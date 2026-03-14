# gnina MCP

> CNN-enhanced molecular docking using Gnina - MCP tools for protein-ligand docking, binding affinity estimation, and virtual screening

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Local Usage (Scripts)](#local-usage-scripts)
- [MCP Server Installation](#mcp-server-installation)
- [Using with Claude Code](#using-with-claude-code)
- [Using with Gemini CLI](#using-with-gemini-cli)
- [Available Tools](#available-tools)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

## Overview

The gnina MCP provides comprehensive molecular docking capabilities using Gnina (CNN-enhanced AutoDock Vina) for protein-ligand interactions. This MCP enables fast scoring, comprehensive docking simulations, virtual screening of compound libraries, and advanced CNN model comparisons - all optimized for drug discovery and cyclic peptide design workflows.

### Features
- **Fast protein-ligand scoring** using empirical and CNN-based scoring functions
- **Comprehensive molecular docking** with pose generation and ranking
- **Virtual screening** for high-throughput compound evaluation
- **CNN model benchmarking** for scoring function optimization
- **Flexible receptor docking** for induced-fit modeling
- **Molecular descriptor analysis** for drug-likeness evaluation
- **Batch processing** for screening multiple ligands against targets

### Directory Structure
```
./
├── README.md               # This file
├── env/                    # Conda environment
├── src/
│   └── server.py           # MCP server
│   └── jobs/               # Job management system
│   └── utils.py            # Shared utilities
├── scripts/
│   ├── score_protein_ligand.py      # Fast protein-ligand scoring
│   ├── dock_ligand.py               # Standard molecular docking
│   ├── virtual_screening.py         # Virtual compound screening
│   ├── compare_cnn_models.py        # CNN model benchmarking
│   ├── flexible_docking.py          # Flexible receptor docking
│   ├── molecular_analysis.py        # Molecular descriptors
│   └── README.md                    # Script documentation
├── examples/
│   └── data/               # Demo data
│       ├── 3rod_rec.pdb    # MDM2 receptor structure
│       ├── 3rod_lig.pdb    # MDM2 ligand
│       ├── 184l_rec.pdb    # Alternative receptor
│       ├── 184l_lig.sdf    # Alternative ligand
│       └── 10gs_*.{pdb,sdf} # Additional test structures
├── configs/                # Configuration files
│   ├── score_config.json   # Scoring parameters
│   ├── dock_config.json    # Docking parameters
│   ├── screen_config.json  # Virtual screening parameters
│   ├── cnn_config.json     # CNN comparison parameters
│   ├── flex_config.json    # Flexible docking parameters
│   └── analysis_config.json # Analysis parameters
└── repo/                   # Original gnina repository
```

---

## Installation

### Prerequisites

**Runtime: conda (chosen deployment runtime)**
- Conda or Mamba (mamba recommended for faster installation)
- Python 3.10+
- Gnina executable (for production use)
- RDKit (installed automatically via conda-forge)

### Quick Setup

```bash
# Navigate to the MCP directory
cd /home/xux/Desktop/AgentMCP/CycPepMCP/tool-mcps/gnina_mcp

# Run the setup script
bash quick_setup.sh
```

### Manual Environment Setup

**Environment Setup (using conda/mamba):**
```bash
cd /home/xux/Desktop/AgentMCP/CycPepMCP/tool-mcps/gnina_mcp

# Create conda environment (use mamba if available)
mamba create -p ./env python=3.10 -y
# or: conda create -p ./env python=3.10 -y

# Activate environment
mamba activate ./env
# or: conda activate ./env

# Install Core Python Packages
mamba run -p ./env pip install loguru click pandas numpy tqdm

# Install FastMCP
mamba run -p ./env pip install --force-reinstall --no-cache-dir fastmcp

# Install RDKit (from conda-forge for better compatibility)
mamba install -p ./env -c conda-forge rdkit -y

# Install Scientific Computing Packages
mamba run -p ./env pip install scipy matplotlib seaborn scikit-learn jupyter pytest
```

**Note:** The chosen runtime for this MCP is **conda**. Use the corresponding instructions above.

---

## Local Usage (Scripts)

You can use the scripts directly without MCP for local processing.

### Available Scripts

| Script | Description | Example |
|--------|-------------|---------|
| `scripts/score_protein_ligand.py` | Score protein-ligand complexes using gnina | See below |
| `scripts/dock_ligand.py` | Molecular docking with CNN scoring | See below |
| `scripts/virtual_screening.py` | High-throughput virtual screening | See below |
| `scripts/compare_cnn_models.py` | CNN model performance comparison | See below |
| `scripts/flexible_docking.py` | Flexible receptor docking | See below |
| `scripts/molecular_analysis.py` | Molecular descriptors and drug-likeness | See below |

### Script Examples

#### Score Protein-Ligand Complex

```bash
# Activate environment
mamba activate ./env

# Run scoring script
python scripts/score_protein_ligand.py \
  --receptor examples/data/3rod_rec.pdb \
  --ligand examples/data/3rod_lig.pdb \
  --output results/scoring_results.csv
```

**Parameters:**
- `--receptor, -r`: Receptor protein structure file (PDB/PDBQT) (required)
- `--ligand, -l`: Ligand structure file (SDF/PDB/PDBQT) (required)
- `--output, -o`: Output CSV file path (optional)
- `--scoring_function`: Gnina scoring function ("default", "vinardo", "ad4_scoring")
- `--cnn_models`: CNN models for scoring (comma-separated)

#### Molecular Docking

```bash
python scripts/dock_ligand.py \
  --receptor examples/data/184l_rec.pdb \
  --ligand examples/data/184l_lig.sdf \
  --autobox_ligand examples/data/184l_lig.sdf \
  --output results/docked_poses.sdf \
  --num_modes 9
```

**Parameters:**
- `--receptor, -r`: Receptor protein structure (PDB) (required)
- `--ligand, -l`: Ligand to dock (SDF/PDB) (required)
- `--autobox_ligand`: Reference ligand for binding site detection
- `--center`: Binding site center (x,y,z coordinates)
- `--size`: Search space size (x,y,z dimensions)
- `--num_modes`: Number of poses to generate (default: 9)
- `--exhaustiveness`: Search exhaustiveness (default: 8)

#### Virtual Screening

```bash
python scripts/virtual_screening.py \
  --receptor examples/data/184l_rec.pdb \
  --ligand_dir examples/data/ \
  --autobox_ligand examples/data/184l_lig.sdf \
  --top_n 50 \
  --output results/screening_results.csv
```

#### Molecular Analysis

```bash
python scripts/molecular_analysis.py \
  --ligand examples/data/184l_lig.sdf \
  --descriptors molecular_weight,logp,tpsa \
  --output results/analysis_results.csv
```

---

## MCP Server Installation

### Option 1: Using fastmcp (Recommended)

```bash
# Install MCP server for Claude Code
fastmcp install src/server.py --name gnina-tools
```

### Option 2: Manual Installation for Claude Code

```bash
# Add MCP server to Claude Code
claude mcp add gnina-tools -- $(pwd)/env/bin/python $(pwd)/src/server.py

# Verify installation
claude mcp list
```

### Option 3: Configure in settings.json

Add to `~/.claude/settings.json`:

```json
{
  "mcpServers": {
    "gnina-tools": {
      "command": "/home/xux/Desktop/AgentMCP/CycPepMCP/tool-mcps/gnina_mcp/env/bin/python",
      "args": ["/home/xux/Desktop/AgentMCP/CycPepMCP/tool-mcps/gnina_mcp/src/server.py"]
    }
  }
}
```

---

## Using with Claude Code

After installing the MCP server, you can use it directly in Claude Code.

### Quick Start

```bash
# Start Claude Code
claude
```

### Example Prompts

#### Tool Discovery
```
What tools are available from gnina-tools?
```

#### Protein-Ligand Scoring (Fast)
```
Score the protein-ligand interaction between @examples/data/3rod_rec.pdb and @examples/data/3rod_lig.pdb using gnina
```

#### Molecular Docking (Submit API)
```
Submit a molecular docking job for receptor @examples/data/184l_rec.pdb and ligand @examples/data/184l_lig.sdf with 9 poses
```

#### Check Job Status
```
Check the status of job abc12345
```

#### Virtual Screening
```
Submit a virtual screening job for receptor @examples/data/184l_rec.pdb against all ligands in @examples/data/ directory, using autobox with @examples/data/184l_lig.sdf as reference, and report top 50 hits
```

#### Molecular Analysis
```
Analyze molecular properties for @examples/data/184l_lig.sdf including molecular weight, LogP, TPSA, and drug-likeness
```

### Using @ References

In Claude Code, use `@` to reference files and directories:

| Reference | Description |
|-----------|-------------|
| `@examples/data/3rod_rec.pdb` | Reference MDM2 receptor structure |
| `@examples/data/184l_lig.sdf` | Reference ligand file |
| `@configs/dock_config.json` | Reference docking configuration |
| `@results/` | Reference output directory |

---

## Using with Gemini CLI

### Configuration

Add to `~/.gemini/settings.json`:

```json
{
  "mcpServers": {
    "gnina-tools": {
      "command": "/home/xux/Desktop/AgentMCP/CycPepMCP/tool-mcps/gnina_mcp/env/bin/python",
      "args": ["/home/xux/Desktop/AgentMCP/CycPepMCP/tool-mcps/gnina_mcp/src/server.py"]
    }
  }
}
```

### Example Prompts

```bash
# Start Gemini CLI
gemini

# Example prompts (same as Claude Code)
> What tools are available from gnina-tools?
> Score protein-ligand complex between receptor.pdb and ligand.sdf
```

---

## Available Tools

### Quick Operations (Sync API)

These tools return results immediately (< 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `score_protein_ligand` | Score protein-ligand complexes using gnina | `receptor_file`, `ligand_file`, `scoring_function`, `cnn_models`, `output_file` |
| `analyze_molecules` | Calculate molecular descriptors and drug-likeness | `ligand_file`, `descriptors`, `output_file` |
| `get_server_info` | Get server capabilities and tool information | None |

### Long-Running Tasks (Submit API)

These tools return a job_id for tracking (> 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `submit_molecular_docking` | Perform molecular docking with pose generation | `receptor_file`, `ligand_file`, `autobox_ligand`, `center`, `size`, `num_modes`, `exhaustiveness`, `cnn_scoring`, `job_name` |
| `submit_virtual_screening` | Screen multiple ligands against target protein | `receptor_file`, `ligand_files`/`ligand_dir`, `autobox_ligand`, `top_n`, `affinity_cutoff`, `job_name` |
| `submit_flexible_docking` | Flexible receptor docking with induced fit | `receptor_file`, `ligand_file`, `flexdist`, `flexdist_ligand`, `flexres`, `compare_rigid`, `job_name` |
| `submit_cnn_comparison` | Benchmark different CNN models | `receptor_file`, `ligand_file`, `models`, `iterations`, `modes`, `job_name` |

### Job Management Tools

| Tool | Description |
|------|-------------|
| `get_job_status` | Check job progress and status |
| `get_job_result` | Get results when job completed |
| `get_job_log` | View job execution logs |
| `cancel_job` | Cancel running job |
| `list_jobs` | List all jobs with filtering |
| `get_queue_info` | Get queue statistics |
| `cleanup_old_jobs` | Clean old completed jobs |

---

## Examples

### Example 1: Quick Protein-Ligand Scoring

**Goal:** Score binding affinity for a protein-ligand complex

**Using Script:**
```bash
python scripts/score_protein_ligand.py \
  --receptor examples/data/3rod_rec.pdb \
  --ligand examples/data/3rod_lig.pdb \
  --output results/mdm2_scoring.csv
```

**Using MCP (in Claude Code):**
```
Score the binding affinity between MDM2 protein @examples/data/3rod_rec.pdb and its ligand @examples/data/3rod_lig.pdb using both default and vinardo scoring functions
```

**Expected Output:**
- Binding affinity scores (kcal/mol) from empirical and CNN models
- Success rate and statistical summary
- Best affinity prediction with confidence

### Example 2: Molecular Docking with Pose Generation

**Goal:** Generate docked poses for drug discovery

**Using Script:**
```bash
python scripts/dock_ligand.py \
  --receptor examples/data/184l_rec.pdb \
  --ligand examples/data/184l_lig.sdf \
  --autobox_ligand examples/data/184l_lig.sdf \
  --output results/docked_poses.sdf \
  --num_modes 9 \
  --exhaustiveness 8
```

**Using MCP (in Claude Code):**
```
Submit molecular docking job:
- Receptor: @examples/data/184l_rec.pdb
- Ligand: @examples/data/184l_lig.sdf
- Use autobox with reference ligand @examples/data/184l_lig.sdf
- Generate 9 poses with exhaustiveness 8
- Enable CNN scoring

Monitor the job and show me the results when complete.
```

### Example 3: Virtual Screening Pipeline

**Goal:** Screen compound library for hit identification

**Using MCP (in Claude Code):**
```
I want to screen all compounds in @examples/data/ directory against receptor @examples/data/184l_rec.pdb:

1. Use @examples/data/184l_lig.sdf as reference for binding site detection
2. Report top 50 hits
3. Use binding affinity cutoff of -7.0 kcal/mol
4. Include CNN-enhanced scoring

Submit the screening job and monitor progress. When complete, analyze the top 10 hits and identify the most promising candidates based on:
- Binding affinity < -8.0 kcal/mol
- CNN score > 0.8
- Drug-like molecular properties
```

### Example 4: CNN Model Benchmarking

**Goal:** Optimize scoring protocol by comparing CNN models

**Using MCP (in Claude Code):**
```
Submit CNN model comparison job for:
- Receptor: @examples/data/3rod_rec.pdb
- Ligand: @examples/data/3rod_lig.pdb
- Compare all available CNN models
- Run 5 iterations for statistical significance
- Test both scoring and docking modes

When complete, recommend the best CNN model for MDM2-ligand interactions based on accuracy and runtime.
```

---

## Demo Data

The `examples/data/` directory contains sample data for testing:

| File | Description | Use With |
|------|-------------|----------|
| `3rod_rec.pdb` | MDM2 receptor structure | All docking and scoring tools |
| `3rod_lig.pdb` | MDM2 small molecule ligand | Scoring and docking |
| `184l_rec.pdb` | Alternative receptor for testing | All tools |
| `184l_lig.sdf` | Alternative ligand in SDF format | All tools |
| `10gs_rec.pdb` | Additional test receptor | Virtual screening |
| `10gs_lig.sdf` | Additional test ligand | Analysis and benchmarking |

---

## Configuration Files

The `configs/` directory contains configuration templates:

| Config | Description | Key Parameters |
|--------|-------------|----------------|
| `dock_config.json` | Docking parameters | `num_modes`, `exhaustiveness`, `cnn_scoring` |
| `score_config.json` | Scoring configuration | `scoring_functions`, `cnn_models`, `timeout` |
| `screen_config.json` | Virtual screening config | `top_n`, `affinity_cutoff`, `workers` |
| `cnn_config.json` | CNN comparison config | `models`, `iterations`, `benchmarking` |
| `flex_config.json` | Flexible docking config | `flexdist`, `flexibility_options` |
| `analysis_config.json` | Molecular analysis config | `descriptors`, `drug_likeness_rules` |

### Config Example

```json
{
  "docking": {
    "num_modes": 9,
    "exhaustiveness": 8,
    "energy_range": 3,
    "spacing": 0.375
  },
  "scoring": {
    "cnn_scoring": "all",
    "cnn_model": "default"
  },
  "binding_site": {
    "size": [20, 20, 20]
  }
}
```

---

## Troubleshooting

### Environment Issues

**Problem:** Environment not found
```bash
# Recreate environment
mamba create -p ./env python=3.10 -y
mamba activate ./env
mamba run -p ./env pip install loguru click pandas numpy tqdm
mamba run -p ./env pip install --force-reinstall --no-cache-dir fastmcp
mamba install -p ./env -c conda-forge rdkit -y
```

**Problem:** RDKit import errors
```bash
# Install RDKit from conda-forge
mamba install -p ./env -c conda-forge rdkit -y
```

**Problem:** Import errors
```bash
# Verify installation
python -c "from src.server import mcp; from rdkit import Chem"
```

### MCP Issues

**Problem:** Server not found in Claude Code
```bash
# Check MCP registration
claude mcp list

# Re-add if needed
claude mcp remove gnina-tools
claude mcp add gnina-tools -- $(pwd)/env/bin/python $(pwd)/src/server.py
```

**Problem:** File not found errors
```
Ensure file paths are absolute or relative to the current working directory.
Use @ references in Claude Code to automatically resolve file paths.
```

**Problem:** Tools not working
```bash
# Test server directly
python -c "
from src.server import mcp
print(list(mcp.list_tools().keys()))
"
```

### Job Issues

**Problem:** Job stuck in pending
```bash
# Check job directory
ls -la jobs/

# View job log
cat jobs/<job_id>/job.log
```

**Problem:** Job failed
```
Use get_job_log with job_id "<job_id>" and tail 100 to see error details
```

**Problem:** Gnina executable not found
```
For testing, a mock gnina executable is provided.
For production, install gnina from: https://github.com/gnina/gnina
```

### Performance Issues

**Problem:** Slow docking performance
```
- Reduce exhaustiveness parameter (default: 8)
- Decrease num_modes (default: 9)
- Use smaller search space size
- Consider using fast CNN models instead of dense models
```

**Problem:** Virtual screening takes too long
```
- Use smaller compound libraries for testing
- Implement more stringent affinity cutoffs
- Consider parallel processing (future enhancement)
```

---

## Development

### Running Tests

```bash
# Activate environment
mamba activate ./env

# Run tests
pytest tests/ -v
```

### Starting Dev Server

```bash
# Run MCP server in dev mode
fastmcp dev src/server.py
```

### Testing Individual Scripts

```bash
# Test scoring script
python scripts/score_protein_ligand.py --help

# Test with demo data
python scripts/score_protein_ligand.py \
  --receptor examples/data/3rod_rec.pdb \
  --ligand examples/data/3rod_lig.pdb
```

---

## License

Based on gnina (Apache 2.0 License)

## Credits

Based on [gnina](https://github.com/gnina/gnina) - CNN-enhanced molecular docking