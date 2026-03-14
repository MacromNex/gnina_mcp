# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

An MCP server wrapping [gnina](https://github.com/gnina/gnina) (CNN-enhanced molecular docking) as tools for protein-ligand scoring, docking, virtual screening, and molecular analysis. Designed for drug discovery and cyclic peptide design workflows.

## Build & Run Commands

```bash
# Setup (download gnina binary into conda env)
bash quick_setup.sh

# Run scripts directly (use env python)
env/bin/python scripts/score_protein_ligand.py -r examples/data/3rod_rec.pdb -l examples/data/3rod_lig.pdb -o results/score.csv
env/bin/python scripts/dock_ligand.py -r examples/data/3rod_rec.pdb -l examples/data/3rod_lig.pdb --autobox_ligand examples/data/3rod_rec_ref.pdb -o results/docked.sdf.gz

# Start MCP server (stdio)
env/bin/python src/server.py

# Docker
docker build -t gnina-mcp .
docker run --gpus all -i gnina-mcp                              # stdio
docker run --gpus all -p 8000:8000 gnina-mcp --transport sse    # SSE

# Tests
env/bin/python tests/run_integration_tests.py
env/bin/python test_mcp_direct.py
```

## Architecture

### Three-Layer Design

```
MCP Server (src/server.py)          ← FastMCP tool definitions
    │
    ├── Sync tools ──→ scripts/*.py  ← Direct function call (score, analyze)
    └── Submit tools ──→ JobManager  ← Background thread + subprocess
                            │
                            └──→ scripts/*.py  (via subprocess.Popen)
```

- **`src/server.py`**: Registers 14 MCP tools via `@mcp.tool()`. Sync tools (fast <10min) import and call script functions directly. Submit tools (long-running) delegate to the job manager.
- **`src/jobs/manager.py`**: Runs scripts as background subprocesses in threads. Each job gets `jobs/{job_id}/` with `metadata.json`, `job.log`, `output.json`. State: PENDING→RUNNING→COMPLETED/FAILED/CANCELLED.
- **`scripts/*.py`**: Self-contained scripts, each with its own `DEFAULT_CONFIG`, CLI parser, and `run_*()` function. Work both standalone and when imported by the MCP server.

### Script ↔ Example Mapping

Each script in `scripts/` is extracted/simplified from a corresponding use case in `examples/`:
- `use_case_1_basic_scoring.py` → `score_protein_ligand.py`
- `use_case_2_standard_docking.py` → `dock_ligand.py`
- `use_case_3_virtual_screening.py` → `virtual_screening.py`
- `use_case_4_cnn_model_comparison.py` → `compare_cnn_models.py`
- `use_case_5_flexible_docking.py` → `flexible_docking.py`
- `use_case_6_python_api.py` → `molecular_analysis.py`

### Config System

Three-level config with deep merge: `DEFAULT_CONFIG` (in script) → config file/dict → CLI/kwargs. All scripts use `_deep_merge()` for recursive dict merging so partial overrides (e.g. `{"gnina": {"verbose": true}}`) don't clobber sibling keys.

### gnina Executable Resolution

Every script and `src/utils.py` implement the same resolution chain:
1. `shutil.which('gnina')` — system PATH
2. `env/bin/gnina` — project conda env (installed by `quick_setup.sh`)
3. `mock_gnina.py` — testing fallback with simulated output
4. `'gnina'` — bare name as last resort

All `subprocess.run()`/`Popen()` calls pass `env=_gnina_subprocess_env()` which prepends `env/lib` to `LD_LIBRARY_PATH` so gnina finds libcudnn.so.9.

## Key Patterns

- **gnina flags**: No `--gpu` flag (gnina auto-detects GPU). Use `--no_gpu` to disable. No `--energy_range`. CNN scoring modes: `none`, `rescore` (default, fast), `refinement` (10x slower), `all` (1000x slower — avoid as default).
- **Config files** in `configs/` are reference templates, not loaded at runtime by default. Each script's `DEFAULT_CONFIG` is authoritative.
- **Three test systems** in `examples/data/`: `3rod` (MDM2), `10gs`, `184l` — each has receptor `.pdb`, ligand `.sdf`/`.pdb`, and reference `.pdb` for autobox.
- **Job manager** uses `fcntl` file locking for atomic metadata writes and `resolve_python_executable()` from `src/utils.py` to find the env python.

## Example Data

| System | Receptor | Ligand | Reference (for autobox) |
|--------|----------|--------|------------------------|
| 3rod | `examples/data/3rod_rec.pdb` | `examples/data/3rod_lig.pdb` | `examples/data/3rod_rec_ref.pdb` |
| 10gs | `examples/data/10gs_rec.pdb` | `examples/data/10gs_lig.sdf` | `examples/data/10gs_rec_ref.pdb` |
| 184l | `examples/data/184l_rec.pdb` | `examples/data/184l_lig.sdf` | `examples/data/184l_rec_ref.pdb` |
