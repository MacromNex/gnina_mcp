# MCP Scripts

Clean, self-contained scripts extracted from use cases for MCP tool wrapping.

## Design Principles

1. **Minimal Dependencies**: Only essential packages imported (pandas, numpy, rdkit)
2. **Self-Contained**: All utility functions inlined
3. **Configurable**: Parameters in config files, not hardcoded
4. **MCP-Ready**: Each script has a main function ready for MCP wrapping

## Scripts

| Script | Description | Dependencies | Config |
|--------|-------------|--------------|--------|
| `score_protein_ligand.py` | Score protein-ligand complexes | pandas | `configs/score_config.json` |
| `dock_ligand.py` | Standard molecular docking | pandas | `configs/dock_config.json` |
| `virtual_screening.py` | High-throughput virtual screening | pandas | `configs/screen_config.json` |
| `compare_cnn_models.py` | CNN model performance comparison | pandas, numpy | `configs/cnn_config.json` |
| `flexible_docking.py` | Flexible receptor docking | pandas | `configs/flex_config.json` |
| `molecular_analysis.py` | Molecular descriptor analysis | pandas, rdkit | `configs/analysis_config.json` |

## Usage

```bash
# Activate environment
mamba activate ./env  # or: conda activate ./env

# Run a script
python scripts/score_protein_ligand.py --input examples/data/3rod_rec.pdb --ligand examples/data/3rod_lig.pdb --output results/output.csv

# With custom config
python scripts/score_protein_ligand.py --config configs/custom.json --input FILE --ligand FILE
```

## For MCP Wrapping (Step 6)

Each script exports a main function that can be wrapped:
```python
from scripts.score_protein_ligand import run_protein_ligand_scoring

# In MCP tool:
@mcp.tool()
def score_protein_ligand_complex(receptor_file: str, ligand_file: str, output_file: str = None):
    return run_protein_ligand_scoring(receptor_file, ligand_file, output_file)
```

## Script Functions

### Core Functions
- `run_protein_ligand_scoring()` - Score existing protein-ligand complexes
- `run_molecular_docking()` - Dock ligand to protein binding site
- `run_virtual_screening()` - Screen multiple ligands against target
- `run_cnn_comparison()` - Compare CNN model performance
- `run_flexible_docking()` - Dock with flexible receptor residues
- `run_molecular_analysis()` - Calculate molecular descriptors and drug-likeness

### Shared Utilities
All scripts include inlined utility functions:
- `parse_gnina_output()` - Parse gnina command line output
- `validate_input_files()` - Validate input file formats
- `run_gnina_command()` - Execute gnina with error handling
- `save_results()` - Save results in various formats