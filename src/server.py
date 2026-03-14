"""MCP Server for Gnina Molecular Docking Tools

Provides both synchronous and asynchronous (submit) APIs for molecular docking
and analysis tools using Gnina (CNN-enhanced molecular docking).
"""

from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import sys

# Setup paths
from utils import setup_paths
paths = setup_paths()
SCRIPT_DIR = paths["script_dir"]
MCP_ROOT = paths["mcp_root"]
SCRIPTS_DIR = paths["scripts_dir"]

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from jobs.manager import job_manager
from loguru import logger

# Create MCP server
mcp = FastMCP("gnina-tools")

# ==============================================================================
# Job Management Tools (for async operations)
# ==============================================================================

@mcp.tool()
def get_job_status(job_id: str) -> dict:
    """
    Get the status of a submitted molecular docking job.

    Args:
        job_id: The job ID returned from a submit_* function

    Returns:
        Dictionary with job status, timestamps, and any errors
    """
    return job_manager.get_job_status(job_id)

@mcp.tool()
def get_job_result(job_id: str) -> dict:
    """
    Get the results of a completed molecular docking job.

    Args:
        job_id: The job ID of a completed job

    Returns:
        Dictionary with the job results or error if not completed
    """
    return job_manager.get_job_result(job_id)

@mcp.tool()
def get_job_log(job_id: str, tail: int = 50) -> dict:
    """
    Get log output from a running or completed job.

    Args:
        job_id: The job ID to get logs for
        tail: Number of lines from end (default: 50, use 0 for all)

    Returns:
        Dictionary with log lines and total line count
    """
    return job_manager.get_job_log(job_id, tail)

@mcp.tool()
def cancel_job(job_id: str) -> dict:
    """
    Cancel a running molecular docking job.

    Args:
        job_id: The job ID to cancel

    Returns:
        Success or error message
    """
    return job_manager.cancel_job(job_id)

@mcp.tool()
def list_jobs(status: Optional[str] = None) -> dict:
    """
    List all submitted molecular docking jobs.

    Args:
        status: Filter by status (pending, running, completed, failed, cancelled)

    Returns:
        List of jobs with their status and metadata
    """
    return job_manager.list_jobs(status)

@mcp.tool()
def get_queue_info() -> dict:
    """
    Get information about the job queue status.

    Returns:
        Dictionary with queue statistics and current state
    """
    return job_manager.get_queue_info()

@mcp.tool()
def cleanup_old_jobs(older_than_days: int = 7) -> dict:
    """
    Clean up completed jobs older than specified days.

    Args:
        older_than_days: Remove completed jobs older than this many days

    Returns:
        Number of jobs cleaned up
    """
    return job_manager.cleanup_old_jobs(older_than_days)

# ==============================================================================
# Synchronous Tools (for fast operations < 10 min)
# ==============================================================================

@mcp.tool()
def score_protein_ligand(
    receptor_file: str,
    ligand_file: str,
    scoring_functions: Optional[List[str]] = None,
    cnn_models: Optional[List[str]] = None,
    output_file: Optional[str] = None
) -> dict:
    """
    Score protein-ligand complexes using gnina scoring functions and CNN models.

    Fast operation - returns results immediately (typically 30 sec - 5 min).

    Args:
        receptor_file: Path to receptor protein structure (PDB/PDBQT)
        ligand_file: Path to ligand structure (SDF/PDB/PDBQT)
        scoring_functions: List of scoring functions (default: ["default", "vinardo"])
        cnn_models: List of CNN models to use for scoring
        output_file: Optional path to save results CSV

    Returns:
        Dictionary with scoring results and statistics
    """
    try:
        from score_protein_ligand import run_protein_ligand_scoring

        config = {}
        if scoring_functions:
            config['scoring'] = {'functions': scoring_functions}
        if cnn_models:
            config.setdefault('scoring', {})['cnn_models'] = cnn_models

        result = run_protein_ligand_scoring(
            receptor_file=receptor_file,
            ligand_file=ligand_file,
            output_file=output_file,
            config=config if config else None
        )
        return {"status": "success", **result}
    except FileNotFoundError as e:
        return {"status": "error", "error": f"File not found: {e}"}
    except ValueError as e:
        return {"status": "error", "error": f"Invalid input: {e}"}
    except Exception as e:
        logger.error(f"Protein-ligand scoring failed: {e}")
        return {"status": "error", "error": str(e)}

@mcp.tool()
def analyze_molecules(
    ligand_file: str,
    output_file: Optional[str] = None
) -> dict:
    """
    Calculate molecular descriptors and drug-likeness for ligand molecules.

    Fast operation using RDKit - returns results immediately (typically 30 sec - 5 min).

    Args:
        ligand_file: Path to ligand structure file (SDF/MOL)
        output_file: Optional path to save analysis results CSV

    Returns:
        Dictionary with molecular descriptors and drug-likeness analysis
    """
    try:
        from molecular_analysis import run_molecular_analysis

        result = run_molecular_analysis(
            ligand_file=ligand_file,
            output_file=output_file
        )
        return {"status": "success", **result}
    except FileNotFoundError as e:
        return {"status": "error", "error": f"File not found: {e}"}
    except ImportError as e:
        return {"status": "error", "error": f"RDKit dependency missing: {e}"}
    except Exception as e:
        logger.error(f"Molecular analysis failed: {e}")
        return {"status": "error", "error": str(e)}

# ==============================================================================
# Submit Tools (for long-running operations > 10 min)
# ==============================================================================

@mcp.tool()
def submit_molecular_docking(
    receptor_file: str,
    ligand_file: str,
    autobox_ligand: Optional[str] = None,
    center: Optional[List[float]] = None,
    size: Optional[List[float]] = None,
    num_modes: int = 9,
    exhaustiveness: int = 8,
    cnn_scoring: bool = True,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit a molecular docking job for background processing.

    This operation typically takes 5-30 minutes. Returns a job_id for tracking.

    Args:
        receptor_file: Path to receptor protein structure (PDB)
        ligand_file: Path to ligand to dock (SDF/PDB)
        autobox_ligand: Path to reference ligand for binding site detection
        center: Binding site center coordinates [x, y, z]
        size: Search space size [x, y, z] (Angstroms)
        num_modes: Number of docking poses to generate
        exhaustiveness: Search exhaustiveness (higher = more thorough)
        cnn_scoring: Use CNN-enhanced scoring
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id. Use:
        - get_job_status(job_id) to check progress
        - get_job_result(job_id) to get results
        - get_job_log(job_id) to see logs
    """
    script_path = str(SCRIPTS_DIR / "dock_ligand.py")

    args = {
        "receptor": receptor_file,
        "ligand": ligand_file,
        "num_modes": num_modes,
        "exhaustiveness": exhaustiveness,
        "cnn_scoring": "rescore" if cnn_scoring else "none"
    }

    # Add optional parameters
    if autobox_ligand:
        args["autobox_ligand"] = autobox_ligand
    if center:
        args["center"] = " ".join(map(str, center))
    if size:
        args["size"] = " ".join(map(str, size))

    return job_manager.submit_job(
        script_path=script_path,
        args=args,
        job_name=job_name or f"docking_{Path(ligand_file).stem}"
    )

@mcp.tool()
def submit_virtual_screening(
    receptor_file: str,
    ligand_files: Optional[List[str]] = None,
    ligand_dir: Optional[str] = None,
    autobox_ligand: Optional[str] = None,
    top_n: int = 100,
    affinity_cutoff: float = -6.0,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit a virtual screening job for background processing.

    Screens multiple ligands against a target protein. This is a long-running
    task (typically 30 minutes to several hours).

    Args:
        receptor_file: Path to target protein structure (PDB)
        ligand_files: List of ligand files to screen
        ligand_dir: Directory containing ligands (alternative to ligand_files)
        autobox_ligand: Reference ligand for binding site detection
        top_n: Number of top hits to report
        affinity_cutoff: Binding affinity cutoff for filtering (kcal/mol)
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id for tracking the screening job
    """
    script_path = str(SCRIPTS_DIR / "virtual_screening.py")

    args = {
        "receptor": receptor_file,
        "top_n": top_n,
        "max_affinity": affinity_cutoff
    }

    # Handle ligand input
    if ligand_files:
        args["ligands"] = " ".join(ligand_files)
    elif ligand_dir:
        args["ligand_dir"] = ligand_dir
    else:
        return {"status": "error", "error": "Must specify either ligand_files or ligand_dir"}

    if autobox_ligand:
        args["autobox_ligand"] = autobox_ligand

    return job_manager.submit_job(
        script_path=script_path,
        args=args,
        job_name=job_name or f"screening_{len(ligand_files) if ligand_files else 'dir'}_ligands"
    )

@mcp.tool()
def submit_flexible_docking(
    receptor_file: str,
    ligand_file: str,
    flexdist: float = 3.5,
    flexdist_ligand: Optional[str] = None,
    flexres: Optional[List[str]] = None,
    compare_rigid: bool = True,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit a flexible receptor docking job for background processing.

    Performs docking with flexible receptor residues. This is a long-running
    task (typically 30-90 minutes).

    Args:
        receptor_file: Path to receptor structure (PDB)
        ligand_file: Path to ligand to dock (SDF/PDB)
        flexdist: Distance cutoff for flexible residues (Angstroms)
        flexdist_ligand: Reference ligand for flexibility selection
        flexres: List of specific residue IDs to make flexible
        compare_rigid: Also perform rigid docking for comparison
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id for tracking the flexible docking job
    """
    script_path = str(SCRIPTS_DIR / "flexible_docking.py")

    args = {
        "receptor": receptor_file,
        "ligand": ligand_file,
        "flexdist": flexdist,
        "compare_rigid": compare_rigid
    }

    if flexdist_ligand:
        args["flexdist_ligand"] = flexdist_ligand
    if flexres:
        args["flexres"] = " ".join(flexres)

    return job_manager.submit_job(
        script_path=script_path,
        args=args,
        job_name=job_name or f"flex_docking_{Path(ligand_file).stem}"
    )

@mcp.tool()
def submit_cnn_comparison(
    receptor_file: str,
    ligand_file: str,
    models: Optional[List[str]] = None,
    iterations: int = 3,
    modes: Optional[List[str]] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit a CNN model comparison job for background processing.

    Benchmarks different CNN models for protein-ligand scoring. This is a
    long-running task (typically 10-60 minutes).

    Args:
        receptor_file: Path to receptor structure (PDB)
        ligand_file: Path to ligand for benchmarking (SDF/PDB)
        models: List of CNN models to compare (default: all available)
        iterations: Number of benchmark iterations per model
        modes: List of scoring modes to test
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id for tracking the comparison job
    """
    script_path = str(SCRIPTS_DIR / "compare_cnn_models.py")

    args = {
        "receptor": receptor_file,
        "ligand": ligand_file,
        "iterations": iterations
    }

    if models:
        args["models"] = " ".join(models)
    if modes:
        args["modes"] = " ".join(modes)

    return job_manager.submit_job(
        script_path=script_path,
        args=args,
        job_name=job_name or f"cnn_comparison_{iterations}_iter"
    )

# ==============================================================================
# Server Information Tools
# ==============================================================================

@mcp.tool()
def get_server_info() -> dict:
    """
    Get information about the Gnina MCP server capabilities.

    Returns:
        Dictionary with server version, available tools, and capabilities
    """
    return {
        "status": "success",
        "server_name": "gnina-tools",
        "version": "1.0.0",
        "description": "CNN-enhanced molecular docking using Gnina",
        "sync_tools": [
            "score_protein_ligand",
            "analyze_molecules"
        ],
        "submit_tools": [
            "submit_molecular_docking",
            "submit_virtual_screening",
            "submit_flexible_docking",
            "submit_cnn_comparison"
        ],
        "job_management": [
            "get_job_status",
            "get_job_result",
            "get_job_log",
            "cancel_job",
            "list_jobs",
            "get_queue_info",
            "cleanup_old_jobs"
        ],
        "supported_formats": {
            "receptor": ["pdb", "pdbqt"],
            "ligand": ["sdf", "pdb", "pdbqt", "mol"],
            "output": ["sdf", "csv", "json"]
        }
    }

# ==============================================================================
# Entry Point
# ==============================================================================

if __name__ == "__main__":
    mcp.run()