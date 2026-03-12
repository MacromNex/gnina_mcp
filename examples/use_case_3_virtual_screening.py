#!/usr/bin/env python3
"""
Use Case 3: Virtual Screening with Gnina
=========================================

This script demonstrates high-throughput virtual screening using gnina for
drug discovery and cyclic peptide optimization. It processes multiple ligands
against a target protein, ranks results by binding affinity and CNN scores,
and provides comprehensive analysis for hit identification.

Use case: Screen large libraries of compounds or cyclic peptides to identify
promising binders for further development.

Requirements:
- gnina executable installed and in PATH
- Receptor PDB file
- Directory containing ligand SDF files or single multi-ligand SDF file
- Optional: Reference ligand for binding site definition

Example Usage:
    python use_case_3_virtual_screening.py --receptor data/184l_rec.pdb --ligand_dir data/ --autobox_ligand data/184l_lig.sdf
    python use_case_3_virtual_screening.py --receptor data/3rod_rec.pdb --ligands ligand1.sdf ligand2.sdf --output screening_results.csv
"""

import argparse
import subprocess
import sys
import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import shutil


def parse_screening_output(output_text):
    """
    Parse gnina output to extract best pose information for screening.

    Args:
        output_text (str): Raw output from gnina

    Returns:
        dict: Best pose information
    """
    output = output_text.decode() if isinstance(output_text, bytes) else output_text

    # Look for the best pose (first line in pose results)
    pose_pattern = r'^\s*1\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+(\S+))?'
    lines = output.split('\n')

    for line in lines:
        match = re.match(pose_pattern, line)
        if match:
            affinity = float(match.group(1))
            rmsd = float(match.group(2)) if match.group(2) != 'inf' else None

            # Handle different output formats
            if len(match.groups()) >= 3:
                cnn_score = float(match.group(3))
            else:
                cnn_score = None

            if len(match.groups()) >= 4 and match.group(4):
                try:
                    cnn_affinity = float(match.group(4))
                except:
                    cnn_affinity = None
            else:
                cnn_affinity = None

            return {
                'affinity': affinity,
                'rmsd': rmsd,
                'cnn_score': cnn_score,
                'cnn_affinity': cnn_affinity
            }

    # If no pose table found, try to extract from summary
    aff_match = re.search(r'Affinity:\s+(\S+)', output)
    cnn_match = re.search(r'CNNscore:\s+(\S+)', output)
    cnn_aff_match = re.search(r'CNNaffinity:\s+(\S+)', output)

    if aff_match:
        return {
            'affinity': float(aff_match.group(1)),
            'rmsd': None,
            'cnn_score': float(cnn_match.group(1)) if cnn_match else None,
            'cnn_affinity': float(cnn_aff_match.group(1)) if cnn_aff_match else None
        }

    return None


def dock_single_ligand(args_tuple):
    """
    Dock a single ligand (for parallel processing).

    Args:
        args_tuple: Tuple of (receptor_path, ligand_path, docking_params, output_dir)

    Returns:
        dict: Docking results for this ligand
    """
    receptor_path, ligand_path, docking_params, output_dir = args_tuple
    ligand_name = Path(ligand_path).stem

    # Create temporary output file if needed
    temp_output = None
    if output_dir:
        temp_output = os.path.join(output_dir, f"{ligand_name}_docked.sdf")

    # Build gnina command
    cmd = ['gnina', '-r', receptor_path, '-l', ligand_path]

    if temp_output:
        cmd.extend(['-o', temp_output])

    # Add docking parameters
    for param, value in docking_params.items():
        if isinstance(value, bool) and value:
            cmd.append(f'--{param}')
        elif not isinstance(value, bool):
            cmd.extend([f'--{param}', str(value)])

    try:
        # Run gnina
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        end_time = time.time()

        output = result.stdout + result.stderr
        pose_info = parse_screening_output(output)

        if pose_info:
            return {
                'ligand_name': ligand_name,
                'ligand_path': ligand_path,
                'affinity': pose_info['affinity'],
                'rmsd': pose_info['rmsd'],
                'cnn_score': pose_info['cnn_score'],
                'cnn_affinity': pose_info['cnn_affinity'],
                'docking_time': end_time - start_time,
                'status': 'success',
                'output_file': temp_output if temp_output and os.path.exists(temp_output) else None
            }
        else:
            return {
                'ligand_name': ligand_name,
                'ligand_path': ligand_path,
                'status': 'failed_parsing',
                'docking_time': end_time - start_time
            }

    except subprocess.CalledProcessError as e:
        return {
            'ligand_name': ligand_name,
            'ligand_path': ligand_path,
            'status': 'failed_docking',
            'error': str(e)
        }
    except Exception as e:
        return {
            'ligand_name': ligand_name,
            'ligand_path': ligand_path,
            'status': 'failed_unknown',
            'error': str(e)
        }


def run_virtual_screening(receptor_path, ligand_paths, docking_params, n_workers=None, output_dir=None):
    """
    Run virtual screening on multiple ligands in parallel.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_paths (list): List of ligand file paths
        docking_params (dict): Gnina docking parameters
        n_workers (int): Number of parallel workers (default: CPU count)
        output_dir (str): Directory to save docked poses (optional)

    Returns:
        pandas.DataFrame: Screening results
    """
    if n_workers is None:
        n_workers = min(mp.cpu_count(), len(ligand_paths))

    print(f"Starting virtual screening with {n_workers} workers...")
    print(f"Screening {len(ligand_paths)} ligands against {Path(receptor_path).name}")

    # Create output directory if needed
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Prepare arguments for parallel processing
    args_list = [(receptor_path, ligand_path, docking_params, output_dir)
                 for ligand_path in ligand_paths]

    results = []
    start_time = time.time()

    # Run docking in parallel
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all tasks
        future_to_ligand = {executor.submit(dock_single_ligand, args): args[1]
                           for args in args_list}

        # Collect results as they complete
        completed = 0
        for future in as_completed(future_to_ligand):
            ligand_path = future_to_ligand[future]
            ligand_name = Path(ligand_path).name

            try:
                result = future.result()
                results.append(result)

                completed += 1
                if result['status'] == 'success':
                    print(f"[{completed}/{len(ligand_paths)}] {ligand_name}: "
                          f"Affinity = {result['affinity']:.3f}, "
                          f"CNN = {result.get('cnn_score', 'N/A'):.3f if result.get('cnn_score') else 'N/A'}")
                else:
                    print(f"[{completed}/{len(ligand_paths)}] {ligand_name}: FAILED ({result['status']})")

            except Exception as e:
                print(f"[{completed}/{len(ligand_paths)}] {ligand_name}: ERROR - {e}")
                results.append({
                    'ligand_name': Path(ligand_path).stem,
                    'ligand_path': ligand_path,
                    'status': 'failed_exception',
                    'error': str(e)
                })

    end_time = time.time()
    total_time = end_time - start_time

    print(f"\nVirtual screening completed in {total_time:.1f} seconds")

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Calculate success rate
    successful = df[df['status'] == 'success']
    print(f"Success rate: {len(successful)}/{len(df)} ({len(successful)/len(df)*100:.1f}%)")

    return df


def analyze_screening_results(results_df):
    """
    Analyze virtual screening results and identify top hits.

    Args:
        results_df (pandas.DataFrame): Screening results

    Returns:
        dict: Analysis summary
    """
    successful = results_df[results_df['status'] == 'success'].copy()

    if len(successful) == 0:
        return {'error': 'No successful docking results'}

    # Basic statistics
    analysis = {
        'total_ligands': len(results_df),
        'successful_ligands': len(successful),
        'success_rate': len(successful) / len(results_df) * 100,
        'mean_docking_time': successful['docking_time'].mean(),
        'total_screening_time': successful['docking_time'].sum()
    }

    # Affinity analysis
    affinities = successful['affinity']
    analysis.update({
        'best_affinity': affinities.min(),
        'worst_affinity': affinities.max(),
        'mean_affinity': affinities.mean(),
        'median_affinity': affinities.median(),
        'affinity_std': affinities.std()
    })

    # CNN score analysis (if available)
    cnn_scores = successful['cnn_score'].dropna()
    if len(cnn_scores) > 0:
        analysis.update({
            'best_cnn_score': cnn_scores.max(),
            'worst_cnn_score': cnn_scores.min(),
            'mean_cnn_score': cnn_scores.mean(),
            'median_cnn_score': cnn_scores.median()
        })

    # Identify top hits
    analysis['top_hits_affinity'] = successful.nsmallest(10, 'affinity')[['ligand_name', 'affinity']].to_dict('records')

    if len(cnn_scores) > 0:
        analysis['top_hits_cnn'] = successful.nlargest(10, 'cnn_score')[['ligand_name', 'cnn_score']].to_dict('records')

    return analysis


def generate_screening_report(results_df, analysis, output_file=None):
    """
    Generate a comprehensive screening report.

    Args:
        results_df (pandas.DataFrame): Screening results
        analysis (dict): Analysis summary
        output_file (str): Output file for report (optional)
    """
    report_lines = []
    report_lines.append("VIRTUAL SCREENING REPORT")
    report_lines.append("=" * 60)

    # Summary statistics
    report_lines.append("\nSCREENING SUMMARY:")
    report_lines.append(f"  Total ligands processed: {analysis['total_ligands']}")
    report_lines.append(f"  Successful dockings: {analysis['successful_ligands']}")
    report_lines.append(f"  Success rate: {analysis['success_rate']:.1f}%")
    report_lines.append(f"  Total screening time: {analysis['total_screening_time']:.1f} seconds")
    report_lines.append(f"  Average docking time: {analysis['mean_docking_time']:.2f} seconds")

    # Affinity statistics
    report_lines.append(f"\nAFFINITY STATISTICS:")
    report_lines.append(f"  Best (most negative): {analysis['best_affinity']:.3f}")
    report_lines.append(f"  Worst (least negative): {analysis['worst_affinity']:.3f}")
    report_lines.append(f"  Mean: {analysis['mean_affinity']:.3f}")
    report_lines.append(f"  Median: {analysis['median_affinity']:.3f}")
    report_lines.append(f"  Standard deviation: {analysis['affinity_std']:.3f}")

    # CNN statistics (if available)
    if 'best_cnn_score' in analysis:
        report_lines.append(f"\nCNN SCORE STATISTICS:")
        report_lines.append(f"  Best (highest): {analysis['best_cnn_score']:.3f}")
        report_lines.append(f"  Worst (lowest): {analysis['worst_cnn_score']:.3f}")
        report_lines.append(f"  Mean: {analysis['mean_cnn_score']:.3f}")
        report_lines.append(f"  Median: {analysis['median_cnn_score']:.3f}")

    # Top hits by affinity
    report_lines.append(f"\nTOP HITS BY AFFINITY:")
    for i, hit in enumerate(analysis['top_hits_affinity'], 1):
        report_lines.append(f"  {i:2d}. {hit['ligand_name']:<20} {hit['affinity']:8.3f}")

    # Top hits by CNN score (if available)
    if 'top_hits_cnn' in analysis:
        report_lines.append(f"\nTOP HITS BY CNN SCORE:")
        for i, hit in enumerate(analysis['top_hits_cnn'], 1):
            report_lines.append(f"  {i:2d}. {hit['ligand_name']:<20} {hit['cnn_score']:8.3f}")

    # Failed ligands summary
    failed = results_df[results_df['status'] != 'success']
    if len(failed) > 0:
        report_lines.append(f"\nFAILED LIGANDS ({len(failed)}):")
        failure_counts = failed['status'].value_counts()
        for status, count in failure_counts.items():
            report_lines.append(f"  {status}: {count}")

    report_text = "\n".join(report_lines)

    # Print to stdout
    print(report_text)

    # Save to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
        print(f"\nReport saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Virtual screening with gnina',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--receptor', '-r', required=True,
                       help='Receptor PDB file')
    parser.add_argument('--ligand_dir',
                       help='Directory containing ligand SDF files')
    parser.add_argument('--ligands', nargs='+',
                       help='List of ligand SDF files')
    parser.add_argument('--output', '-o',
                       help='Output CSV file for results')
    parser.add_argument('--output_dir',
                       help='Directory to save docked poses')
    parser.add_argument('--autobox_ligand',
                       help='Reference ligand for automatic box definition')
    parser.add_argument('--center', nargs=3, type=float,
                       help='Box center coordinates (x y z)')
    parser.add_argument('--size', nargs=3, type=float, default=[20, 20, 20],
                       help='Box size (x y z) in Angstroms (default: 20 20 20)')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                       help='Exhaustiveness of search (default: 8)')
    parser.add_argument('--num_modes', type=int, default=1,
                       help='Number of binding modes per ligand (default: 1)')
    parser.add_argument('--cnn_scoring',
                       choices=['none', 'rescore', 'refinement'],
                       default='rescore',
                       help='CNN scoring mode (default: rescore)')
    parser.add_argument('--scoring',
                       choices=['default', 'vinardo', 'ad4_scoring'],
                       help='Empirical scoring function')
    parser.add_argument('--workers', '-w', type=int,
                       help='Number of parallel workers (default: CPU count)')
    parser.add_argument('--top_n', type=int, default=10,
                       help='Number of top hits to highlight (default: 10)')
    parser.add_argument('--no_gpu', action='store_true',
                       help='Disable GPU acceleration')
    parser.add_argument('--report',
                       help='Save detailed report to file')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.receptor):
        print(f"Error: Receptor file {args.receptor} not found")
        sys.exit(1)

    # Collect ligand files
    ligand_paths = []
    if args.ligand_dir:
        ligand_dir = Path(args.ligand_dir)
        if not ligand_dir.exists():
            print(f"Error: Ligand directory {args.ligand_dir} not found")
            sys.exit(1)
        ligand_paths = list(ligand_dir.glob('*.sdf'))
    elif args.ligands:
        ligand_paths = [Path(f) for f in args.ligands if Path(f).exists()]
        missing = [f for f in args.ligands if not Path(f).exists()]
        if missing:
            print(f"Error: Missing ligand files: {missing}")
            sys.exit(1)
    else:
        print("Error: Must specify either --ligand_dir or --ligands")
        sys.exit(1)

    if not ligand_paths:
        print("Error: No valid ligand files found")
        sys.exit(1)

    # Setup docking parameters
    docking_params = {
        'exhaustiveness': args.exhaustiveness,
        'num_modes': args.num_modes,
        'cnn_scoring': args.cnn_scoring,
        'seed': 42  # For reproducibility
    }

    if args.autobox_ligand:
        if not os.path.exists(args.autobox_ligand):
            print(f"Error: Autobox ligand file {args.autobox_ligand} not found")
            sys.exit(1)
        docking_params['autobox_ligand'] = args.autobox_ligand
    elif args.center:
        docking_params.update({
            'center_x': args.center[0],
            'center_y': args.center[1],
            'center_z': args.center[2],
            'size_x': args.size[0],
            'size_y': args.size[1],
            'size_z': args.size[2]
        })
    else:
        print("Error: Must specify either --autobox_ligand or --center")
        sys.exit(1)

    if args.scoring:
        docking_params['scoring'] = args.scoring

    if args.no_gpu:
        docking_params['no_gpu'] = True

    print(f"Starting virtual screening:")
    print(f"  Receptor: {args.receptor}")
    print(f"  Ligands: {len(ligand_paths)} files")
    print(f"  Workers: {args.workers or mp.cpu_count()}")
    print(f"  CNN scoring: {args.cnn_scoring}")

    # Run virtual screening
    results_df = run_virtual_screening(
        args.receptor,
        [str(p) for p in ligand_paths],
        docking_params,
        args.workers,
        args.output_dir
    )

    # Analyze results
    analysis = analyze_screening_results(results_df)

    if 'error' not in analysis:
        # Generate report
        generate_screening_report(results_df, analysis, args.report)

        # Save results to CSV
        if args.output:
            # Sort by affinity for CSV output
            successful = results_df[results_df['status'] == 'success'].copy()
            if len(successful) > 0:
                successful_sorted = successful.sort_values('affinity')
                successful_sorted.to_csv(args.output, index=False)
                print(f"\nDetailed results saved to: {args.output}")

        print(f"\nVirtual screening completed successfully!")
        print(f"Top hit: {analysis['top_hits_affinity'][0]['ligand_name']} "
              f"(affinity: {analysis['top_hits_affinity'][0]['affinity']:.3f})")

    else:
        print(f"Virtual screening failed: {analysis['error']}")
        sys.exit(1)


if __name__ == '__main__':
    main()