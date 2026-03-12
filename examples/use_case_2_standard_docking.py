#!/usr/bin/env python3
"""
Use Case 2: Standard Molecular Docking with Gnina
==================================================

This script demonstrates how to perform molecular docking using gnina with
CNN-enhanced scoring. It includes automatic binding site detection, multiple
pose generation, and various CNN scoring modes.

Use case: Dock small molecules or cyclic peptides into protein binding sites
with deep learning-enhanced scoring for improved pose prediction and ranking.

Requirements:
- gnina executable installed and in PATH
- Receptor PDB file
- Ligand SDF/PDB file
- Optional: Reference ligand for autobox definition

Example Usage:
    python use_case_2_standard_docking.py --receptor data/184l_rec.pdb --ligand data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf
    python use_case_2_standard_docking.py --receptor data/3rod_rec.pdb --ligand data/3rod_lig.pdb --output docked_poses.sdf --num_modes 20
"""

import argparse
import subprocess
import sys
import os
import re
import pandas as pd
from pathlib import Path
import time


def parse_docking_output(output_text):
    """
    Parse gnina docking output to extract pose information.

    Args:
        output_text (str): Raw output from gnina

    Returns:
        list: List of dictionaries containing pose information
    """
    output = output_text.decode() if isinstance(output_text, bytes) else output_text
    poses = []

    # Look for pose table in output (mode, affinity, rmsd, cnn_score, cnn_affinity)
    pose_pattern = r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$'
    lines = output.split('\n')

    in_pose_section = False
    for line in lines:
        # Check if we're in the pose results section
        if 'mode |' in line or 'affinity' in line.lower():
            in_pose_section = True
            continue

        if in_pose_section:
            match = re.match(pose_pattern, line)
            if match:
                mode = int(match.group(1))
                affinity = float(match.group(2))

                # Handle different output formats (sometimes has RMSD, sometimes not)
                try:
                    if len(match.groups()) == 4:
                        rmsd = float(match.group(3))
                        cnn_score = float(match.group(4))
                        cnn_affinity = None
                    elif len(match.groups()) == 5:
                        rmsd = float(match.group(3))
                        cnn_score = float(match.group(4))
                        cnn_affinity = float(match.group(5))
                    else:
                        rmsd = None
                        cnn_score = float(match.group(3))
                        cnn_affinity = None
                except:
                    rmsd = None
                    cnn_score = None
                    cnn_affinity = None

                poses.append({
                    'mode': mode,
                    'affinity': affinity,
                    'rmsd': rmsd,
                    'cnn_score': cnn_score,
                    'cnn_affinity': cnn_affinity
                })
            elif line.strip() == '':
                # End of pose section
                break

    return poses


def run_gnina_docking(receptor_path, ligand_path, output_path=None, **kwargs):
    """
    Run gnina molecular docking.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_path (str): Path to ligand SDF/PDB file
        output_path (str): Output file for docked poses (optional)
        **kwargs: Additional gnina parameters

    Returns:
        tuple: (poses_info, raw_output, success)
    """
    # Build gnina command
    cmd = ['gnina', '-r', receptor_path, '-l', ligand_path]

    # Add output file if specified
    if output_path:
        cmd.extend(['-o', output_path])

    # Add optional parameters
    param_mapping = {
        'autobox_ligand': '--autobox_ligand',
        'autobox_add': '--autobox_add',
        'center_x': '--center_x',
        'center_y': '--center_y',
        'center_z': '--center_z',
        'size_x': '--size_x',
        'size_y': '--size_y',
        'size_z': '--size_z',
        'num_modes': '--num_modes',
        'exhaustiveness': '--exhaustiveness',
        'seed': '--seed',
        'cnn': '--cnn',
        'cnn_scoring': '--cnn_scoring',
        'scoring': '--scoring',
        'pose_sort_order': '--pose_sort_order',
        'min_rmsd_filter': '--min_rmsd_filter',
        'no_gpu': '--no_gpu'
    }

    for param, flag in param_mapping.items():
        if param in kwargs and kwargs[param] is not None:
            if isinstance(kwargs[param], bool) and kwargs[param]:
                cmd.append(flag)
            elif not isinstance(kwargs[param], bool):
                cmd.extend([flag, str(kwargs[param])])

    try:
        # Run gnina
        print(f"Running command: {' '.join(cmd)}")
        start_time = time.time()

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        end_time = time.time()
        runtime = end_time - start_time

        output = result.stdout + result.stderr

        # Parse pose information
        poses = parse_docking_output(output)

        print(f"Docking completed in {runtime:.1f} seconds")
        print(f"Generated {len(poses)} poses")

        return poses, output, True

    except subprocess.CalledProcessError as e:
        print(f"Error running gnina: {e}")
        print(f"Command: {' '.join(cmd)}")
        print(f"Error output: {e.stderr}")
        return [], None, False
    except FileNotFoundError:
        print("Error: gnina executable not found. Please ensure gnina is installed and in PATH.")
        return [], None, False


def compare_docking_modes(receptor_path, ligand_path, autobox_ligand=None):
    """
    Compare different CNN scoring modes during docking.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_path (str): Path to ligand file
        autobox_ligand (str): Path to reference ligand for autoboxing

    Returns:
        dict: Results for each scoring mode
    """
    scoring_modes = {
        'rescore': 'rescore',
        'refinement': 'refinement',
        'none_vinardo': 'none'
    }

    results = {}

    for mode_name, cnn_scoring in scoring_modes.items():
        print(f"\n{'='*50}")
        print(f"Testing {mode_name.upper()} mode")
        print(f"{'='*50}")

        kwargs = {
            'num_modes': 9,
            'exhaustiveness': 8,
            'seed': 42,
            'cnn_scoring': cnn_scoring
        }

        if autobox_ligand:
            kwargs['autobox_ligand'] = autobox_ligand

        if mode_name == 'none_vinardo':
            kwargs['scoring'] = 'vinardo'

        poses, output, success = run_gnina_docking(
            receptor_path, ligand_path, **kwargs
        )

        results[mode_name] = {
            'poses': poses,
            'success': success,
            'output': output
        }

        if success and poses:
            best_pose = min(poses, key=lambda x: x['affinity'])
            print(f"Best pose - Affinity: {best_pose['affinity']:.3f}, CNN Score: {best_pose.get('cnn_score', 'N/A')}")

    return results


def analyze_pose_quality(poses, reference_pose=None):
    """
    Analyze the quality and diversity of generated poses.

    Args:
        poses (list): List of pose dictionaries
        reference_pose (dict): Reference pose for comparison (optional)

    Returns:
        dict: Analysis results
    """
    if not poses:
        return {}

    affinities = [p['affinity'] for p in poses]
    cnn_scores = [p['cnn_score'] for p in poses if p['cnn_score'] is not None]
    rmsds = [p['rmsd'] for p in poses if p['rmsd'] is not None]

    analysis = {
        'num_poses': len(poses),
        'affinity_range': max(affinities) - min(affinities),
        'best_affinity': min(affinities),
        'worst_affinity': max(affinities),
        'mean_affinity': sum(affinities) / len(affinities)
    }

    if cnn_scores:
        analysis.update({
            'best_cnn_score': max(cnn_scores),
            'worst_cnn_score': min(cnn_scores),
            'mean_cnn_score': sum(cnn_scores) / len(cnn_scores)
        })

    if rmsds:
        analysis.update({
            'mean_rmsd': sum(rmsds) / len(rmsds),
            'max_rmsd': max(rmsds),
            'min_rmsd': min(rmsds)
        })

    return analysis


def main():
    parser = argparse.ArgumentParser(
        description='Perform molecular docking using gnina',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--receptor', '-r', required=True,
                       help='Receptor PDB file')
    parser.add_argument('--ligand', '-l', required=True,
                       help='Ligand SDF/PDB file')
    parser.add_argument('--output', '-o',
                       help='Output SDF file for docked poses')
    parser.add_argument('--autobox_ligand',
                       help='Reference ligand for automatic box definition')
    parser.add_argument('--center', nargs=3, type=float,
                       help='Box center coordinates (x y z)')
    parser.add_argument('--size', nargs=3, type=float, default=[20, 20, 20],
                       help='Box size (x y z) in Angstroms (default: 20 20 20)')
    parser.add_argument('--num_modes', type=int, default=9,
                       help='Number of binding modes to generate (default: 9)')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                       help='Exhaustiveness of search (default: 8)')
    parser.add_argument('--cnn_scoring',
                       choices=['none', 'rescore', 'refinement', 'all'],
                       default='rescore',
                       help='CNN scoring mode (default: rescore)')
    parser.add_argument('--cnn_model',
                       choices=['default', 'fast', 'dense', 'default1.0', 'general_default2018'],
                       default='default',
                       help='CNN model to use (default: default)')
    parser.add_argument('--scoring',
                       choices=['default', 'vinardo', 'ad4_scoring', 'dkoes_fast'],
                       help='Empirical scoring function')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--no_gpu', action='store_true',
                       help='Disable GPU acceleration')
    parser.add_argument('--compare_modes', action='store_true',
                       help='Compare different CNN scoring modes')
    parser.add_argument('--pose_sort_order',
                       choices=['CNNscore', 'CNNaffinity', 'Energy'],
                       help='How to sort output poses')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.receptor):
        print(f"Error: Receptor file {args.receptor} not found")
        sys.exit(1)

    if not os.path.exists(args.ligand):
        print(f"Error: Ligand file {args.ligand} not found")
        sys.exit(1)

    if args.autobox_ligand and not os.path.exists(args.autobox_ligand):
        print(f"Error: Autobox ligand file {args.autobox_ligand} not found")
        sys.exit(1)

    print(f"Receptor: {args.receptor}")
    print(f"Ligand: {args.ligand}")
    print(f"Output: {args.output or 'stdout only'}")

    if args.compare_modes:
        print("\nComparing different CNN scoring modes...")
        results = compare_docking_modes(args.receptor, args.ligand, args.autobox_ligand)

        print("\n" + "="*60)
        print("COMPARISON SUMMARY")
        print("="*60)

        for mode_name, result in results.items():
            if result['success']:
                analysis = analyze_pose_quality(result['poses'])
                print(f"\n{mode_name.upper()}:")
                print(f"  Poses generated: {analysis['num_poses']}")
                print(f"  Best affinity: {analysis['best_affinity']:.3f}")
                print(f"  Affinity range: {analysis['affinity_range']:.3f}")
                if 'best_cnn_score' in analysis:
                    print(f"  Best CNN score: {analysis['best_cnn_score']:.3f}")
            else:
                print(f"\n{mode_name.upper()}: FAILED")

    else:
        # Standard docking
        kwargs = {
            'num_modes': args.num_modes,
            'exhaustiveness': args.exhaustiveness,
            'seed': args.seed,
            'cnn_scoring': args.cnn_scoring,
            'no_gpu': args.no_gpu
        }

        if args.autobox_ligand:
            kwargs['autobox_ligand'] = args.autobox_ligand
        elif args.center:
            kwargs['center_x'] = args.center[0]
            kwargs['center_y'] = args.center[1]
            kwargs['center_z'] = args.center[2]
            kwargs['size_x'] = args.size[0]
            kwargs['size_y'] = args.size[1]
            kwargs['size_z'] = args.size[2]
        else:
            print("Error: Must specify either --autobox_ligand or --center for binding site definition")
            sys.exit(1)

        if args.cnn_model != 'default':
            kwargs['cnn'] = args.cnn_model

        if args.scoring:
            kwargs['scoring'] = args.scoring

        if args.pose_sort_order:
            kwargs['pose_sort_order'] = args.pose_sort_order

        poses, output, success = run_gnina_docking(
            args.receptor, args.ligand, args.output, **kwargs
        )

        if success:
            analysis = analyze_pose_quality(poses)

            print("\n" + "="*60)
            print("DOCKING RESULTS")
            print("="*60)

            if poses:
                # Display pose table
                print(f"\n{'Mode':<6}{'Affinity':<10}{'RMSD':<8}{'CNN Score':<12}{'CNN Affinity':<12}")
                print("-" * 50)

                for pose in poses:
                    rmsd_str = f"{pose['rmsd']:.3f}" if pose['rmsd'] is not None else "N/A"
                    cnn_score_str = f"{pose['cnn_score']:.3f}" if pose['cnn_score'] is not None else "N/A"
                    cnn_aff_str = f"{pose['cnn_affinity']:.3f}" if pose['cnn_affinity'] is not None else "N/A"

                    print(f"{pose['mode']:<6}{pose['affinity']:<10.3f}{rmsd_str:<8}{cnn_score_str:<12}{cnn_aff_str:<12}")

                print(f"\nAnalysis Summary:")
                print(f"  Total poses: {analysis['num_poses']}")
                print(f"  Best affinity: {analysis['best_affinity']:.3f}")
                print(f"  Worst affinity: {analysis['worst_affinity']:.3f}")
                print(f"  Mean affinity: {analysis['mean_affinity']:.3f}")

                if 'best_cnn_score' in analysis:
                    print(f"  Best CNN score: {analysis['best_cnn_score']:.3f}")
                    print(f"  Mean CNN score: {analysis['mean_cnn_score']:.3f}")

                if args.output:
                    print(f"\nDocked poses saved to: {args.output}")

            else:
                print("No poses were generated successfully.")
        else:
            print("Docking failed. Please check your input files and parameters.")
            sys.exit(1)


if __name__ == '__main__':
    main()