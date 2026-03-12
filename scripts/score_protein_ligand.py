#!/usr/bin/env python3
"""
Script: score_protein_ligand.py
Description: Score protein-ligand complexes using gnina scoring functions and CNN models

Original Use Case: examples/use_case_1_basic_scoring.py
Dependencies Removed: None (already minimal)

Usage:
    python scripts/score_protein_ligand.py --receptor <receptor.pdb> --ligand <ligand.sdf> --output <output.csv>

Example:
    python scripts/score_protein_ligand.py --receptor examples/data/3rod_rec.pdb --ligand examples/data/3rod_lig.pdb --output results/output.csv
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
from pathlib import Path
from typing import Union, Optional, Dict, Any, List, Tuple

# Essential scientific packages
import pandas as pd

# ==============================================================================
# Configuration (extracted from use case)
# ==============================================================================
DEFAULT_CONFIG = {
    "scoring": {
        "functions": ["default", "vinardo"],
        "cnn_models": ["default"],
        "score_only": True,
        "batch_size": 10
    },
    "output": {
        "format": "csv",
        "include_raw_output": False,
        "float_precision": 3
    },
    "gnina": {
        "executable": "gnina",
        "timeout": 300,
        "verbose": False
    }
}

# ==============================================================================
# Inlined Utility Functions (simplified from repo)
# ==============================================================================
def parse_gnina_output(output_text: str) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Parse gnina output to extract affinity and CNN scores. Inlined from original use case."""
    output = output_text.decode() if isinstance(output_text, bytes) else output_text

    affinity = None
    cnn_score = None
    cnn_affinity = None

    # Extract Vina affinity score
    match = re.search(r'Affinity:\s+(\S+)', output)
    if match:
        try:
            affinity = float(match.group(1))
        except ValueError:
            pass

    # Extract CNN score
    match = re.search(r'CNNscore:\s+(\S+)', output)
    if match:
        try:
            cnn_score = float(match.group(1))
        except ValueError:
            pass

    # Extract CNN affinity
    match = re.search(r'CNNaffinity:\s+(\S+)', output)
    if match:
        try:
            cnn_affinity = float(match.group(1))
        except ValueError:
            pass

    return affinity, cnn_score, cnn_affinity


def validate_input_files(receptor_path: Union[str, Path], ligand_path: Union[str, Path]) -> bool:
    """Validate that input files exist and have correct extensions."""
    receptor_path = Path(receptor_path)
    ligand_path = Path(ligand_path)

    if not receptor_path.exists():
        raise FileNotFoundError(f"Receptor file not found: {receptor_path}")

    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_path}")

    # Check file extensions
    valid_receptor_exts = {'.pdb', '.pdbqt'}
    valid_ligand_exts = {'.sdf', '.pdb', '.pdbqt', '.mol2'}

    if receptor_path.suffix.lower() not in valid_receptor_exts:
        raise ValueError(f"Invalid receptor file format. Expected: {valid_receptor_exts}")

    if ligand_path.suffix.lower() not in valid_ligand_exts:
        raise ValueError(f"Invalid ligand file format. Expected: {valid_ligand_exts}")

    return True


def run_gnina_command(receptor_path: str, ligand_path: str, scoring_function: str = 'default',
                     cnn_model: str = 'default', gnina_executable: str = 'gnina',
                     timeout: int = 300) -> Tuple[Optional[float], Optional[float], Optional[float], str]:
    """Execute gnina scoring command with error handling."""
    # Build gnina command
    cmd = [gnina_executable, '-r', str(receptor_path), '-l', str(ligand_path), '--score_only']

    # Add scoring function if not default
    if scoring_function != 'default':
        cmd.extend(['--scoring', scoring_function])

    # Add CNN model if not default
    if cnn_model != 'default':
        cmd.extend(['--cnn', cnn_model])

    try:
        # Run gnina
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        output = result.stdout + result.stderr

        # Parse scores
        affinity, cnn_score, cnn_affinity = parse_gnina_output(output)

        return affinity, cnn_score, cnn_affinity, output

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Gnina command timed out after {timeout} seconds")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Gnina command failed: {e.stderr}")
    except FileNotFoundError:
        raise RuntimeError(f"Gnina executable not found: {gnina_executable}")


def save_results(results_data: List[Dict], output_path: Union[str, Path], format_type: str = "csv") -> None:
    """Save results in specified format."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format_type.lower() == "csv":
        df = pd.DataFrame(results_data)
        df.to_csv(output_path, index=False, float_format='%.3f')
    elif format_type.lower() == "json":
        with open(output_path, 'w') as f:
            json.dump(results_data, f, indent=2)
    else:
        raise ValueError(f"Unsupported output format: {format_type}")


# ==============================================================================
# Core Function (main logic extracted from use case)
# ==============================================================================
def run_protein_ligand_scoring(
    receptor_file: Union[str, Path],
    ligand_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Score protein-ligand complex using gnina scoring functions and CNN models.

    Args:
        receptor_file: Path to receptor PDB file
        ligand_file: Path to ligand file (SDF, PDB, etc.)
        output_file: Path to save results (optional)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - results: List of scoring results
            - output_file: Path to output file (if saved)
            - metadata: Execution metadata

    Example:
        >>> result = run_protein_ligand_scoring("receptor.pdb", "ligand.sdf", "output.csv")
        >>> print(result['results'][0]['affinity'])
    """
    # Setup
    receptor_file = Path(receptor_file)
    ligand_file = Path(ligand_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Validate inputs
    validate_input_files(receptor_file, ligand_file)

    # Get scoring parameters
    scoring_functions = config['scoring']['functions']
    cnn_models = config['scoring']['cnn_models']
    gnina_executable = config['gnina']['executable']
    timeout = config['gnina']['timeout']

    # Run scoring for all function/model combinations
    results = []
    for scoring_func in scoring_functions:
        for cnn_model in cnn_models:
            try:
                affinity, cnn_score, cnn_affinity, raw_output = run_gnina_command(
                    receptor_file, ligand_file, scoring_func, cnn_model,
                    gnina_executable, timeout
                )

                result_entry = {
                    "receptor": str(receptor_file.name),
                    "ligand": str(ligand_file.name),
                    "scoring_function": scoring_func,
                    "cnn_model": cnn_model,
                    "affinity": affinity,
                    "cnn_score": cnn_score,
                    "cnn_affinity": cnn_affinity,
                    "success": True
                }

                if config['output']['include_raw_output']:
                    result_entry["raw_output"] = raw_output

                results.append(result_entry)

            except Exception as e:
                # Record failed attempts
                results.append({
                    "receptor": str(receptor_file.name),
                    "ligand": str(ligand_file.name),
                    "scoring_function": scoring_func,
                    "cnn_model": cnn_model,
                    "affinity": None,
                    "cnn_score": None,
                    "cnn_affinity": None,
                    "success": False,
                    "error": str(e)
                })

    # Save output if requested
    output_path = None
    if output_file:
        output_path = Path(output_file)
        save_results(results, output_path, config['output']['format'])

    # Calculate summary statistics
    successful_results = [r for r in results if r['success']]
    stats = {}
    if successful_results:
        affinities = [r['affinity'] for r in successful_results if r['affinity'] is not None]
        if affinities:
            stats = {
                "best_affinity": min(affinities),
                "worst_affinity": max(affinities),
                "mean_affinity": sum(affinities) / len(affinities),
                "successful_runs": len(successful_results),
                "total_runs": len(results)
            }

    return {
        "results": results,
        "output_file": str(output_path) if output_path else None,
        "statistics": stats,
        "metadata": {
            "receptor_file": str(receptor_file),
            "ligand_file": str(ligand_file),
            "config": config,
            "total_combinations": len(scoring_functions) * len(cnn_models),
            "successful": len(successful_results),
            "failed": len(results) - len(successful_results)
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
    parser.add_argument('--ligand', '-l', required=True, help='Ligand file (SDF, PDB, etc.)')
    parser.add_argument('--output', '-o', help='Output CSV file path')
    parser.add_argument('--config', '-c', help='Config file (JSON)')
    parser.add_argument('--scoring_functions', nargs='+',
                       default=['default'],
                       choices=['default', 'vinardo', 'ad4_scoring', 'dkoes_fast'],
                       help='Scoring functions to test')
    parser.add_argument('--cnn_models', nargs='+',
                       default=['default'],
                       choices=['default', 'fast', 'dense', 'default1.0', 'general_default2018'],
                       help='CNN models to test')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    # Override config with command line arguments
    config_overrides = {}
    if args.scoring_functions != ['default']:
        config_overrides['scoring'] = config_overrides.get('scoring', {})
        config_overrides['scoring']['functions'] = args.scoring_functions
    if args.cnn_models != ['default']:
        config_overrides['scoring'] = config_overrides.get('scoring', {})
        config_overrides['scoring']['cnn_models'] = args.cnn_models
    if args.verbose:
        config_overrides['gnina'] = config_overrides.get('gnina', {})
        config_overrides['gnina']['verbose'] = True

    # Merge configs
    if config:
        config.update(config_overrides)
    else:
        config = config_overrides

    try:
        # Run scoring
        result = run_protein_ligand_scoring(
            receptor_file=args.receptor,
            ligand_file=args.ligand,
            output_file=args.output,
            config=config
        )

        # Print results
        print("="*60)
        print("PROTEIN-LIGAND SCORING RESULTS")
        print("="*60)

        if result['results']:
            df = pd.DataFrame(result['results'])
            print(df.to_string(index=False, float_format='%.3f'))

        if result['statistics']:
            print("\n" + "="*60)
            print("SUMMARY STATISTICS")
            print("="*60)
            stats = result['statistics']
            print(f"Best affinity: {stats.get('best_affinity', 'N/A'):.3f}")
            print(f"Mean affinity: {stats.get('mean_affinity', 'N/A'):.3f}")
            print(f"Success rate: {stats.get('successful_runs', 0)}/{stats.get('total_runs', 0)}")

        if result['output_file']:
            print(f"\nResults saved to: {result['output_file']}")

        return result

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()