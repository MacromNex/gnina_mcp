# gnina MCP Server

**CNN-enhanced molecular docking, scoring, and virtual screening via Docker**

An MCP (Model Context Protocol) server for gnina molecular docking with 14 tools:
- Score protein-ligand complexes with empirical and CNN scoring functions
- Dock ligands into protein binding sites with pose generation
- Screen compound libraries via high-throughput virtual screening
- Compare CNN model performance for scoring optimization
- Perform flexible receptor docking for induced-fit modeling
- Analyze molecular descriptors and drug-likeness

## Quick Start with Docker

### Approach 1: Pull Pre-built Image from GitHub

The fastest way to get started. A pre-built Docker image is automatically published to GitHub Container Registry on every release.

```bash
# Pull the latest image
docker pull ghcr.io/macromnex/gnina_mcp:latest

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add gnina-tools -- docker run -i --rm --user `id -u`:`id -g` --gpus all -v `pwd`:`pwd` ghcr.io/macromnex/gnina_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker with GPU support (`nvidia-docker` or Docker with NVIDIA runtime)
- Claude Code installed

That's it! The gnina MCP server is now available in Claude Code.

---

### Approach 2: Build Docker Image Locally

Build the image yourself and install it into Claude Code. Useful for customization or offline environments.

```bash
# Clone the repository
git clone https://github.com/MacromNex/gnina_mcp.git
cd gnina_mcp

# Build the Docker image
docker build -t gnina_mcp:latest .

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add gnina-tools -- docker run -i --rm --user `id -u`:`id -g` --gpus all -v `pwd`:`pwd` gnina_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker with GPU support
- Claude Code installed
- Git (to clone the repository)

**About the Docker Flags:**
- `-i` — Interactive mode for Claude Code
- `--rm` — Automatically remove container after exit
- `` --user `id -u`:`id -g` `` — Runs the container as your current user, so output files are owned by you (not root)
- `--gpus all` — Grants access to all available GPUs
- `-v` — Mounts your project directory so the container can access your data

---

## Verify Installation

After adding the MCP server, you can verify it's working:

```bash
# List registered MCP servers
claude mcp list

# You should see 'gnina-tools' in the output
```

In Claude Code, you can now use all 14 gnina tools:

**Sync tools (fast, return immediately):**
- `score_protein_ligand` — Score protein-ligand binding affinity
- `analyze_molecules` — Calculate molecular descriptors and drug-likeness
- `get_server_info` — Get server capabilities

**Submit tools (long-running, return job_id):**
- `submit_molecular_docking` — Dock ligands with pose generation
- `submit_virtual_screening` — Screen compound libraries
- `submit_flexible_docking` — Flexible receptor docking
- `submit_cnn_comparison` — Benchmark CNN scoring models

**Job management:**
- `get_job_status`, `get_job_result`, `get_job_log`, `cancel_job`, `list_jobs`, `get_queue_info`, `cleanup_old_jobs`

---

## Usage Examples

Once registered, you can use the gnina tools directly in Claude Code. Here are some common workflows:

### Example 1: Protein-Ligand Scoring

```
I have a protein receptor at /path/to/receptor.pdb and a ligand at /path/to/ligand.sdf. Can you score their binding affinity using score_protein_ligand with both default and vinardo scoring functions?
```

### Example 2: Molecular Docking

```
Dock the ligand /path/to/ligand.sdf into the binding site of /path/to/receptor.pdb, using /path/to/ref_ligand.pdb as the reference for autobox. Generate 9 poses and save results to /path/to/docked.sdf.gz.
```

### Example 3: Virtual Screening

```
Screen all SDF files in /path/to/ligands/ against receptor /path/to/receptor.pdb. Use /path/to/ref.pdb for binding site detection. Report the top 50 hits with affinity better than -7.0 kcal/mol.
```

### Example 4: CNN Model Comparison

```
Compare the default and fast CNN models for scoring the complex of /path/to/receptor.pdb and /path/to/ligand.sdf. Run 3 iterations and recommend the best model based on accuracy and speed.
```

---

## Next Steps

- **Detailed documentation**: See [detail.md](detail.md) for comprehensive guides on:
  - Available MCP tools and parameters
  - Local Python environment setup (alternative to Docker)
  - Script usage and CLI options
  - Configuration file formats
  - Example workflows and troubleshooting

---

## Troubleshooting

**Docker not found?**
```bash
docker --version  # Install Docker if missing
```

**GPU not accessible?**
- Ensure NVIDIA Docker runtime is installed
- Check with `docker run --gpus all ubuntu nvidia-smi`
- CPU mode works without GPU (gnina uses `--no_gpu` automatically)

**Claude Code not found?**
```bash
# Install Claude Code
npm install -g @anthropic-ai/claude-code
```

---

## License

Based on [gnina](https://github.com/gnina/gnina) (dual licensed under GPL and Apache)
