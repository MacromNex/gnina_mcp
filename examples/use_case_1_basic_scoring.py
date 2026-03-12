#!/usr/bin/env python3
"""
Use Case 1: Basic Scoring with Gnina
====================================

This script demonstrates how to use gnina for scoring protein-ligand complexes
without performing docking. It evaluates binding affinity using both empirical
scoring functions and CNN-based scoring.

Use case: Quickly score existing protein-ligand poses for virtual screening
or binding affinity estimation.

Requirements:
- gnina executable installed and in PATH
- Receptor PDB file
- Ligand SDF file (positioned in binding site)

Example Usage:
    python use_case_1_basic_scoring.py --receptor data/3rod_rec.pdb --ligand data/3rod_lig.pdb
    python use_case_1_basic_scoring.py --receptor data/184l_rec.pdb --ligand data/184l_lig.sdf --output scoring_results.csv
"""

import argparse
import subprocess
import sys
import os
import re
import pandas as pd
from pathlib import Path


def parse_gnina_output(output_text):
    """
    Parse gnina output to extract affinity and CNN scores.

    Args:
        output_text (str): Raw output from gnina

    Returns:
        tuple: (affinity, cnn_score, cnn_affinity)
    """
    output = output_text.decode() if isinstance(output_text, bytes) else output_text

    affinity = None
    cnn_score = None
    cnn_affinity = None

    # Extract Vina affinity score
    match = re.search(r'Affinity:\s+(\S+)', output)
    if match:
        affinity = float(match.group(1))

    # Extract CNN score
    match = re.search(r'CNNscore:\s+(\S+)', output)
    if match:
        cnn_score = float(match.group(1))

    # Extract CNN affinity
    match = re.search(r'CNNaffinity:\s+(\S+)', output)
    if match:
        cnn_affinity = float(match.group(1))

    return affinity, cnn_score, cnn_affinity


def run_gnina_scoring(receptor_path, ligand_path, scoring_function='default', cnn_model='default'):
    """
    Run gnina in score-only mode.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_path (str): Path to ligand SDF/PDB file
        scoring_function (str): Empirical scoring function ('default', 'vinardo', 'ad4_scoring')
        cnn_model (str): CNN model to use ('default', 'fast', 'dense', 'default1.0')

    Returns:
        tuple: (affinity, cnn_score, cnn_affinity, raw_output)
    """
    # Build gnina command
    cmd = ['gnina', '-r', receptor_path, '-l', ligand_path, '--score_only']

    # Add scoring function if not default
    if scoring_function != 'default':
        cmd.extend(['--scoring', scoring_function])

    # Add CNN model if not default
    if cnn_model != 'default':
        cmd.extend(['--cnn', cnn_model])

    try:
        # Run gnina
        print(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stdout + result.stderr

        # Parse scores
        affinity, cnn_score, cnn_affinity = parse_gnina_output(output)

        return affinity, cnn_score, cnn_affinity, output

    except subprocess.CalledProcessError as e:
        print(f"Error running gnina: {e}")
        print(f"Command: {' '.join(cmd)}")
        print(f"Error output: {e.stderr}")
        return None, None, None, None
    except FileNotFoundError:
        print("Error: gnina executable not found. Please ensure gnina is installed and in PATH.")
        return None, None, None, None


def score_multiple_complexes(receptor_path, ligand_paths, scoring_functions=['default'], cnn_models=['default']):
    """
    Score multiple ligands against a receptor with different scoring functions.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_paths (list): List of ligand file paths
        scoring_functions (list): List of scoring functions to test
        cnn_models (list): List of CNN models to test

    Returns:
        pandas.DataFrame: Results table with all scores
    """
    results = []

    for ligand_path in ligand_paths:
        ligand_name = Path(ligand_path).stem

        for scoring_func in scoring_functions:
            for cnn_model in cnn_models:
                print(f"\nScoring {ligand_name} with {scoring_func} scoring and {cnn_model} CNN...")

                affinity, cnn_score, cnn_affinity, raw_output = run_gnina_scoring(
                    receptor_path, ligand_path, scoring_func, cnn_model
                )

                if affinity is not None:
                    results.append({
                        'ligand': ligand_name,
                        'receptor': Path(receptor_path).stem,
                        'scoring_function': scoring_func,
                        'cnn_model': cnn_model,
                        'affinity': affinity,
                        'cnn_score': cnn_score,
                        'cnn_affinity': cnn_affinity
                    })

                    print(f"  Affinity: {affinity:.3f}")
                    print(f"  CNN Score: {cnn_score:.3f}")
                    print(f"  CNN Affinity: {cnn_affinity:.3f}")
                else:
                    print(f"  Failed to score {ligand_name}")

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description='Score protein-ligand complexes using gnina',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--receptor', '-r', required=True,
                       help='Receptor PDB file')
    parser.add_argument('--ligand', '-l',
                       help='Ligand SDF/PDB file (for single ligand scoring)')
    parser.add_argument('--ligand_dir',
                       help='Directory containing multiple ligand files')
    parser.add_argument('--output', '-o',
                       help='Output CSV file for results (optional)')
    parser.add_argument('--scoring_functions', nargs='+',
                       default=['default', 'vinardo'],
                       choices=['default', 'vinardo', 'ad4_scoring', 'dkoes_fast'],
                       help='Empirical scoring functions to test')
    parser.add_argument('--cnn_models', nargs='+',
                       default=['default'],
                       choices=['default', 'fast', 'dense', 'default1.0', 'general_default2018'],
                       help='CNN models to test')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Print detailed output')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.receptor):
        print(f"Error: Receptor file {args.receptor} not found")
        sys.exit(1)

    # Determine ligand files
    ligand_paths = []
    if args.ligand:
        if not os.path.exists(args.ligand):
            print(f"Error: Ligand file {args.ligand} not found")
            sys.exit(1)
        ligand_paths = [args.ligand]
    elif args.ligand_dir:
        ligand_dir = Path(args.ligand_dir)
        if not ligand_dir.exists():
            print(f"Error: Ligand directory {args.ligand_dir} not found")
            sys.exit(1)
        ligand_paths = list(ligand_dir.glob('*.sdf')) + list(ligand_dir.glob('*.pdb'))
        if not ligand_paths:
            print(f"Error: No ligand files found in {args.ligand_dir}")
            sys.exit(1)
    else:
        print("Error: Must specify either --ligand or --ligand_dir")
        sys.exit(1)

    print(f"Receptor: {args.receptor}")
    print(f"Ligands: {len(ligand_paths)} files")
    print(f"Scoring functions: {args.scoring_functions}")
    print(f"CNN models: {args.cnn_models}")

    # Run scoring
    results_df = score_multiple_complexes(
        args.receptor,
        [str(p) for p in ligand_paths],
        args.scoring_functions,
        args.cnn_models
    )

    if len(results_df) > 0:
        print("\n" + "="*60)
        print("SCORING RESULTS")
        print("="*60)
        print(results_df.to_string(index=False, float_format='%.3f'))

        # Save results
        if args.output:
            results_df.to_csv(args.output, index=False)
            print(f"\nResults saved to: {args.output}")

        # Summary statistics
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)

        for scoring_func in args.scoring_functions:
            subset = results_df[results_df['scoring_function'] == scoring_func]
            if len(subset) > 0:
                print(f"\n{scoring_func.upper()} SCORING:")
                print(f"  Best affinity: {subset['affinity'].min():.3f}")
                print(f"  Worst affinity: {subset['affinity'].max():.3f}")
                print(f"  Mean affinity: {subset['affinity'].mean():.3f}")

                if subset['cnn_affinity'].notna().any():
                    print(f"  Best CNN affinity: {subset['cnn_affinity'].max():.3f}")
                    print(f"  Worst CNN affinity: {subset['cnn_affinity'].min():.3f}")
                    print(f"  Mean CNN affinity: {subset['cnn_affinity'].mean():.3f}")

    else:
        print("No successful scoring results obtained.")
        sys.exit(1)


if __name__ == '__main__':
    main()