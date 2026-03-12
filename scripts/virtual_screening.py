#!/usr/bin/env python3
"""
Script: virtual_screening.py
Description: High-throughput virtual screening using gnina for drug discovery

Original Use Case: examples/use_case_3_virtual_screening.py
Dependencies Removed: multiprocessing (simplified to sequential processing for MCP compatibility)

Usage:
    python scripts/virtual_screening.py --receptor <receptor.pdb> --ligand_dir <ligands_dir> --output <results.csv>

Example:
    python scripts/virtual_screening.py --receptor examples/data/184l_rec.pdb --ligands examples/data/*.sdf --output results/screening.csv
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import subprocess
import sys
import os
import re
import json
import glob
from pathlib import Path
from typing import Union, Optional, Dict, Any, List, Tuple

# Essential scientific packages
import pandas as pd

# ==============================================================================
# Configuration (extracted from use case)
# ==============================================================================
DEFAULT_CONFIG = {
    "screening": {
        "workers": 1,  # Simplified to 1 for MCP compatibility
        "top_n": 10,
        "batch_size": 100,
        "continue_on_error": True
    },
    "docking": {
        "num_modes": 3,
        "exhaustiveness": 4,
        "energy_range": 2,
        "autobox_add": 4
    },
    "filtering": {
        "max_affinity": 0,
        "min_cnn_score": -1000,
        "remove_duplicates": True
    },
    "output": {
        "format": "csv",
        "save_poses": False,
        "include_failed": True,
        "detailed_report": True
    },
    "gnina": {
        "executable": "gnina",
        "timeout": 600,
        "gpu": True,
        "verbose": False
    }
}

# ==============================================================================
# Inlined Utility Functions
# ==============================================================================
def parse_screening_output(output_text: str, ligand_name: str) -> Dict[str, Any]:
    """Parse gnina output for virtual screening results."""
    output = output_text.decode() if isinstance(output_text, bytes) else output_text

    result = {
        "ligand": ligand_name,
        "success": False,
        "affinity": None,
        "cnn_score": None,
        "cnn_affinity": None,
        "poses_generated": 0,
        "error": None
    }

    try:
        # Extract best affinity (first pose)
        affinity_match = re.search(r'1\s+(-?\d+\.?\d*)', output)
        if affinity_match:
            result["affinity"] = float(affinity_match.group(1))
            result["success"] = True

        # Extract CNN scores
        cnn_score_match = re.search(r'CNNscore:\s+(\S+)', output)
        if cnn_score_match:
            try:
                result["cnn_score"] = float(cnn_score_match.group(1))
            except ValueError:
                pass

        cnn_affinity_match = re.search(r'CNNaffinity:\s+(\S+)', output)
        if cnn_affinity_match:
            try:
                result["cnn_affinity"] = float(cnn_affinity_match.group(1))
            except ValueError:
                pass

        # Count poses generated
        pose_matches = re.findall(r'^\s*\d+\s+', output, re.MULTILINE)
        result["poses_generated"] = len(pose_matches)

    except Exception as e:
        result["error"] = str(e)

    return result


def validate_screening_inputs(receptor_path: Union[str, Path],
                            ligand_paths: List[Union[str, Path]],
                            autobox_ligand: Optional[Union[str, Path]] = None) -> bool:
    """Validate virtual screening input parameters."""
    receptor_path = Path(receptor_path)

    if not receptor_path.exists():
        raise FileNotFoundError(f"Receptor file not found: {receptor_path}")

    if not ligand_paths:
        raise ValueError("No ligand files provided")

    # Check ligand files exist
    valid_ligands = []
    for ligand_path in ligand_paths:
        ligand_path = Path(ligand_path)
        if ligand_path.exists():
            valid_ligands.append(ligand_path)
        else:
            print(f"Warning: Ligand file not found: {ligand_path}")

    if not valid_ligands:
        raise ValueError("No valid ligand files found")

    # Check autobox ligand
    if autobox_ligand:
        autobox_path = Path(autobox_ligand)
        if not autobox_path.exists():
            raise FileNotFoundError(f"Autobox ligand file not found: {autobox_path}")

    return True


def screen_single_ligand(receptor_path: str, ligand_path: str, autobox_ligand: str,
                        config: Dict[str, Any]) -> Dict[str, Any]:
    """Screen a single ligand against the receptor."""
    ligand_path = Path(ligand_path)

    # Build gnina command
    cmd = [
        config['gnina']['executable'],
        '-r', str(receptor_path),
        '-l', str(ligand_path),
        '--num_modes', str(config['docking']['num_modes']),
        '--exhaustiveness', str(config['docking']['exhaustiveness']),
        '--energy_range', str(config['docking']['energy_range'])
    ]

    # Add autobox if specified
    if autobox_ligand:
        cmd.extend(['--autobox_ligand', str(autobox_ligand)])
        cmd.extend(['--autobox_add', str(config['docking']['autobox_add'])])

    # Add CNN scoring
    cmd.extend(['--cnn_scoring', 'all'])

    # Add GPU if available
    if config['gnina']['gpu']:
        cmd.append('--gpu')

    try:
        # Run gnina
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True,
            timeout=config['gnina']['timeout']
        )
        output = result.stdout + result.stderr

        # Parse results
        return parse_screening_output(output, ligand_path.name)

    except subprocess.TimeoutExpired:
        return {
            "ligand": ligand_path.name,
            "success": False,
            "error": f"Timeout after {config['gnina']['timeout']} seconds",
            "affinity": None,
            "cnn_score": None,
            "cnn_affinity": None,
            "poses_generated": 0
        }
    except subprocess.CalledProcessError as e:
        return {
            "ligand": ligand_path.name,
            "success": False,
            "error": f"Gnina error: {e.stderr}",
            "affinity": None,
            "cnn_score": None,
            "cnn_affinity": None,
            "poses_generated": 0
        }
    except Exception as e:
        return {
            "ligand": ligand_path.name,
            "success": False,
            "error": str(e),
            "affinity": None,
            "cnn_score": None,
            "cnn_affinity": None,
            "poses_generated": 0
        }


def filter_and_rank_results(results: List[Dict[str, Any]], config: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Filter and rank screening results."""
    # Filter successful results
    successful_results = [r for r in results if r['success'] and r['affinity'] is not None]

    # Apply filtering criteria
    filtered_results = []
    for result in successful_results:
        # Filter by affinity
        if result['affinity'] > config['filtering']['max_affinity']:
            continue

        # Filter by CNN score if available
        if (result.get('cnn_score') is not None and
            result['cnn_score'] < config['filtering']['min_cnn_score']):
            continue

        filtered_results.append(result)

    # Sort by affinity (best first)
    filtered_results.sort(key=lambda x: x['affinity'])

    # Return top N results
    top_n = config['screening']['top_n']
    return filtered_results[:top_n]


def generate_screening_report(results: List[Dict[str, Any]],
                            filtered_results: List[Dict[str, Any]]) -> str:
    """Generate a comprehensive screening report."""
    total_ligands = len(results)
    successful = len([r for r in results if r['success']])
    failed = total_ligands - successful

    report = []
    report.append("="*60)
    report.append("VIRTUAL SCREENING REPORT")
    report.append("="*60)
    report.append(f"Total ligands screened: {total_ligands}")
    report.append(f"Successful: {successful}")
    report.append(f"Failed: {failed}")
    report.append(f"Success rate: {(successful/total_ligands)*100:.1f}%")
    report.append("")

    if successful > 0:
        affinities = [r['affinity'] for r in results if r['success'] and r['affinity'] is not None]
        report.append("AFFINITY STATISTICS:")
        report.append(f"  Best affinity: {min(affinities):.3f} kcal/mol")
        report.append(f"  Worst affinity: {max(affinities):.3f} kcal/mol")
        report.append(f"  Mean affinity: {sum(affinities)/len(affinities):.3f} kcal/mol")
        report.append("")

    report.append(f"TOP {len(filtered_results)} HITS:")
    for i, hit in enumerate(filtered_results, 1):
        report.append(f"  {i}. {hit['ligand']}: {hit['affinity']:.3f} kcal/mol")
        if hit.get('cnn_affinity'):
            report.append(f"     CNN affinity: {hit['cnn_affinity']:.3f}")

    return "\n".join(report)


# ==============================================================================
# Core Function
# ==============================================================================
def run_virtual_screening(
    receptor_file: Union[str, Path],
    ligand_files: List[Union[str, Path]] = None,
    ligand_dir: Optional[Union[str, Path]] = None,
    output_file: Optional[Union[str, Path]] = None,
    autobox_ligand: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform high-throughput virtual screening using gnina.

    Args:
        receptor_file: Path to receptor PDB file
        ligand_files: List of ligand files to screen
        ligand_dir: Directory containing ligand files (alternative to ligand_files)
        output_file: Path to save results CSV
        autobox_ligand: Reference ligand for binding site detection
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - results: All screening results
            - top_hits: Filtered and ranked top hits
            - report: Screening report text
            - output_file: Path to output file
            - metadata: Execution metadata

    Example:
        >>> result = run_virtual_screening("receptor.pdb",
        ...                               ligand_files=["lig1.sdf", "lig2.sdf"],
        ...                               autobox_ligand="ref.sdf",
        ...                               output_file="screening.csv")
        >>> print(f"Found {len(result['top_hits'])} hits")
    """
    # Setup
    receptor_file = Path(receptor_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Get ligand files
    ligand_paths = []
    if ligand_files:
        ligand_paths = [Path(f) for f in ligand_files]
    elif ligand_dir:
        ligand_dir = Path(ligand_dir)
        # Find all ligand files in directory
        for ext in ['*.sdf', '*.pdb', '*.mol2', '*.pdbqt']:
            ligand_paths.extend(ligand_dir.glob(ext))
    else:
        raise ValueError("Either ligand_files or ligand_dir must be provided")

    # Validate inputs
    validate_screening_inputs(receptor_file, ligand_paths, autobox_ligand)

    print(f"Starting virtual screening of {len(ligand_paths)} ligands...")

    # Screen ligands sequentially (simplified for MCP)
    results = []
    for i, ligand_path in enumerate(ligand_paths, 1):
        if config['gnina']['verbose']:
            print(f"Screening {i}/{len(ligand_paths)}: {ligand_path.name}")

        result = screen_single_ligand(
            receptor_file, ligand_path, autobox_ligand, config
        )
        results.append(result)

        # Stop early on critical errors if configured
        if not config['screening']['continue_on_error'] and not result['success']:
            print(f"Stopping due to error: {result['error']}")
            break

    # Filter and rank results
    top_hits = filter_and_rank_results(results, config)

    # Generate report
    report = generate_screening_report(results, top_hits)

    # Save results
    output_path = None
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save detailed results
        df = pd.DataFrame(results)
        df.to_csv(output_path, index=False, float_format='%.3f')

        # Save top hits separately
        if top_hits:
            hits_path = output_path.parent / f"{output_path.stem}_top_hits.csv"
            hits_df = pd.DataFrame(top_hits)
            hits_df.to_csv(hits_path, index=False, float_format='%.3f')

        # Save report
        if config['output']['detailed_report']:
            report_path = output_path.parent / f"{output_path.stem}_report.txt"
            with open(report_path, 'w') as f:
                f.write(report)

    return {
        "results": results,
        "top_hits": top_hits,
        "report": report,
        "output_file": str(output_path) if output_path else None,
        "metadata": {
            "receptor_file": str(receptor_file),
            "total_ligands": len(ligand_paths),
            "successful_screens": len([r for r in results if r['success']]),
            "config": config
        }
    }


# ==============================================================================
# CLI Interface
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--receptor', '-r', required=True, help='Receptor PDB file')
    parser.add_argument('--ligand_dir', help='Directory containing ligand files')
    parser.add_argument('--ligands', nargs='*', help='Specific ligand files to screen')
    parser.add_argument('--output', '-o', help='Output CSV file')
    parser.add_argument('--autobox_ligand', help='Reference ligand for autobox')
    parser.add_argument('--top_n', type=int, default=10, help='Number of top hits to report')
    parser.add_argument('--max_affinity', type=float, default=0,
                       help='Maximum affinity cutoff (kcal/mol)')
    parser.add_argument('--timeout', type=int, default=600,
                       help='Timeout per ligand (seconds)')
    parser.add_argument('--config', '-c', help='Config file (JSON)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Validate ligand input
    ligand_files = []
    if args.ligands:
        ligand_files = args.ligands
    elif args.ligand_dir:
        # Expand wildcards in ligand_dir
        ligand_dir = Path(args.ligand_dir)
        for ext in ['*.sdf', '*.pdb', '*.mol2']:
            ligand_files.extend(glob.glob(str(ligand_dir / ext)))
    else:
        print("Error: Must specify either --ligands or --ligand_dir")
        sys.exit(1)

    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    # Override config with command line arguments
    config_overrides = {
        "screening": {
            "top_n": args.top_n
        },
        "filtering": {
            "max_affinity": args.max_affinity
        },
        "gnina": {
            "timeout": args.timeout,
            "verbose": args.verbose
        }
    }

    # Merge configs
    if config:
        for section, params in config_overrides.items():
            if section not in config:
                config[section] = {}
            config[section].update(params)
    else:
        config = config_overrides

    try:
        # Run virtual screening
        result = run_virtual_screening(
            receptor_file=args.receptor,
            ligand_files=ligand_files,
            output_file=args.output,
            autobox_ligand=args.autobox_ligand,
            config=config
        )

        # Print results
        print(result['report'])

        if result['output_file']:
            print(f"\nDetailed results saved to: {result['output_file']}")

        return result

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()