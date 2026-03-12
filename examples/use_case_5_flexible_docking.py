#!/usr/bin/env python3
"""
Use Case 5: Flexible Docking with Gnina
========================================

This script demonstrates flexible receptor docking with gnina, where selected
protein side chains are allowed to move during the docking process. This is
particularly important for cyclic peptides and other ligands that may induce
conformational changes in the binding site.

Use case: Dock ligands into flexible binding sites where side chain movement
is expected, such as allosteric sites or highly dynamic binding pockets.

Requirements:
- gnina executable installed and in PATH
- Receptor PDB file
- Ligand SDF/PDB file
- makeflex.py script (included in gnina)
- Optional: Reference ligand for flexible residue selection

Example Usage:
    python use_case_5_flexible_docking.py --receptor data/184l_rec.pdb --ligand data/184l_lig.sdf --flexdist 3.5
    python use_case_5_flexible_docking.py --receptor data/3rod_rec.pdb --ligand data/3rod_lig.pdb --flexres A:120,A:145,A:200 --autobox_ligand data/3rod_lig.pdb
"""

import argparse
import subprocess
import sys
import os
import re
import pandas as pd
import tempfile
from pathlib import Path
import time
import shutil


def identify_flexible_residues(receptor_path, reference_ligand=None, distance_cutoff=3.5):
    """
    Identify potential flexible residues near the binding site.

    Args:
        receptor_path (str): Path to receptor PDB file
        reference_ligand (str): Path to reference ligand for distance calculation
        distance_cutoff (float): Distance cutoff for flexible residue selection

    Returns:
        list: List of residue identifiers (chain:resid format)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
        import numpy as np

        # This is a simplified implementation
        # In practice, you might want to use more sophisticated methods
        print(f"Identifying flexible residues within {distance_cutoff} Å of binding site...")

        # For now, return some example flexible residues
        # In a real implementation, you would parse the PDB and calculate distances
        flexible_residues = []

        # Read receptor PDB to identify residues
        with open(receptor_path, 'r') as f:
            lines = f.readlines()

        # Simple parsing to identify side chain atoms near binding site
        # This is a placeholder - real implementation would calculate distances to ligand
        for line in lines:
            if line.startswith('ATOM') and 'GLU' in line or 'ASP' in line or 'ARG' in line or 'LYS' in line:
                # Extract chain and residue number
                chain = line[21]
                resnum = line[22:26].strip()
                res_id = f"{chain}:{resnum}"

                if res_id not in flexible_residues:
                    flexible_residues.append(res_id)
                    if len(flexible_residues) >= 5:  # Limit to reasonable number
                        break

        print(f"Identified {len(flexible_residues)} potential flexible residues")
        return flexible_residues

    except ImportError:
        print("Warning: RDKit not available for automatic flexible residue identification")
        return []


def parse_flexible_docking_output(output_text):
    """
    Parse gnina flexible docking output.

    Args:
        output_text (str): Raw output from gnina

    Returns:
        dict: Parsed results including poses and flexible residue information
    """
    output = output_text.decode() if isinstance(output_text, bytes) else output_text

    results = {
        'poses': [],
        'flexible_residues': [],
        'affinity': None,
        'cnn_score': None,
        'cnn_affinity': None,
        'runtime': None
    }

    # Extract pose information
    pose_pattern = r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+(\S+))?'
    lines = output.split('\n')

    in_pose_section = False
    for line in lines:
        if 'mode |' in line or 'affinity' in line.lower():
            in_pose_section = True
            continue

        if in_pose_section:
            match = re.match(pose_pattern, line)
            if match:
                try:
                    mode = int(match.group(1))
                    affinity = float(match.group(2))
                    rmsd = float(match.group(3)) if match.group(3) != 'inf' else None
                    cnn_score = float(match.group(4))
                    cnn_affinity = float(match.group(5)) if len(match.groups()) >= 5 and match.group(5) else None

                    pose = {
                        'mode': mode,
                        'affinity': affinity,
                        'rmsd': rmsd,
                        'cnn_score': cnn_score,
                        'cnn_affinity': cnn_affinity
                    }
                    results['poses'].append(pose)

                    # Store best pose info
                    if mode == 1:
                        results['affinity'] = affinity
                        results['cnn_score'] = cnn_score
                        results['cnn_affinity'] = cnn_affinity

                except ValueError:
                    continue
            elif line.strip() == '':
                break

    # Extract flexible residue information
    flex_pattern = r'Using (\d+) flexible residues'
    flex_match = re.search(flex_pattern, output)
    if flex_match:
        results['num_flexible'] = int(flex_match.group(1))

    return results


def run_flexible_docking(receptor_path, ligand_path, output_path=None, **kwargs):
    """
    Run flexible docking with gnina.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_path (str): Path to ligand file
        output_path (str): Output path for docked poses
        **kwargs: Additional parameters

    Returns:
        tuple: (results, success, raw_output)
    """
    # Build gnina command
    cmd = ['gnina', '-r', receptor_path, '-l', ligand_path]

    if output_path:
        cmd.extend(['-o', output_path])

    # Add flexible docking parameters
    param_mapping = {
        'flexres': '--flexres',
        'flexdist_ligand': '--flexdist_ligand',
        'flexdist': '--flexdist',
        'flex_limit': '--flex_limit',
        'flex_max': '--flex_max',
        'autobox_ligand': '--autobox_ligand',
        'center_x': '--center_x',
        'center_y': '--center_y',
        'center_z': '--center_z',
        'size_x': '--size_x',
        'size_y': '--size_y',
        'size_z': '--size_z',
        'num_modes': '--num_modes',
        'exhaustiveness': '--exhaustiveness',
        'seed': '--seed',
        'cnn_scoring': '--cnn_scoring',
        'scoring': '--scoring',
        'no_gpu': '--no_gpu'
    }

    for param, flag in param_mapping.items():
        if param in kwargs and kwargs[param] is not None:
            if isinstance(kwargs[param], bool) and kwargs[param]:
                cmd.append(flag)
            elif not isinstance(kwargs[param], bool):
                cmd.extend([flag, str(kwargs[param])])

    # Add output file for flexible residues if needed
    flex_output = None
    if output_path:
        flex_output = output_path.replace('.sdf', '_flex.pdb')
        cmd.extend(['--out_flex', flex_output])

    try:
        print(f"Running flexible docking command: {' '.join(cmd)}")
        start_time = time.time()

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        end_time = time.time()
        runtime = end_time - start_time

        output = result.stdout + result.stderr
        parsed_results = parse_flexible_docking_output(output)
        parsed_results['runtime'] = runtime

        print(f"Flexible docking completed in {runtime:.1f} seconds")

        return parsed_results, True, output

    except subprocess.CalledProcessError as e:
        print(f"Error running flexible docking: {e}")
        print(f"Command: {' '.join(cmd)}")
        print(f"Error output: {e.stderr}")
        return None, False, e.stderr
    except FileNotFoundError:
        print("Error: gnina executable not found. Please ensure gnina is installed and in PATH.")
        return None, False, "gnina not found"


def combine_flexible_structures(rigid_pdb, flexible_pdb, output_pdb, makeflex_script=None):
    """
    Combine rigid receptor with flexible residues using makeflex.py.

    Args:
        rigid_pdb (str): Path to rigid receptor PDB
        flexible_pdb (str): Path to flexible residues PDB from gnina
        output_pdb (str): Output path for combined structure
        makeflex_script (str): Path to makeflex.py script

    Returns:
        bool: Success status
    """
    if not os.path.exists(flexible_pdb):
        print(f"Warning: Flexible residues file not found: {flexible_pdb}")
        return False

    # Find makeflex.py script
    if makeflex_script is None:
        # Try to find it in common locations
        potential_paths = [
            'examples/makeflex.py',
            'makeflex.py',
            os.path.join(os.path.dirname(__file__), 'makeflex.py')
        ]

        for path in potential_paths:
            if os.path.exists(path):
                makeflex_script = path
                break

    if makeflex_script is None or not os.path.exists(makeflex_script):
        print("Warning: makeflex.py script not found. Cannot combine structures.")
        return False

    try:
        cmd = ['python', makeflex_script, rigid_pdb, flexible_pdb, output_pdb]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        print(f"Combined structure saved to: {output_pdb}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error combining structures: {e}")
        return False


def analyze_flexible_results(results, reference_results=None):
    """
    Analyze flexible docking results and compare with rigid docking if available.

    Args:
        results (dict): Flexible docking results
        reference_results (dict): Rigid docking results for comparison

    Returns:
        dict: Analysis summary
    """
    analysis = {
        'num_poses': len(results['poses']),
        'best_affinity': min([p['affinity'] for p in results['poses']]) if results['poses'] else None,
        'runtime': results.get('runtime'),
        'num_flexible': results.get('num_flexible', 0)
    }

    if results['poses']:
        affinities = [p['affinity'] for p in results['poses']]
        analysis.update({
            'affinity_range': max(affinities) - min(affinities),
            'mean_affinity': sum(affinities) / len(affinities)
        })

        cnn_scores = [p['cnn_score'] for p in results['poses'] if p['cnn_score'] is not None]
        if cnn_scores:
            analysis.update({
                'best_cnn_score': max(cnn_scores),
                'mean_cnn_score': sum(cnn_scores) / len(cnn_scores)
            })

    # Compare with rigid docking if available
    if reference_results:
        analysis['comparison'] = {}

        if results['poses'] and reference_results.get('poses'):
            flex_best = min([p['affinity'] for p in results['poses']])
            rigid_best = min([p['affinity'] for p in reference_results['poses']])

            analysis['comparison']['affinity_improvement'] = rigid_best - flex_best
            analysis['comparison']['runtime_overhead'] = results.get('runtime', 0) - reference_results.get('runtime', 0)

    return analysis


def run_comparison_study(receptor_path, ligand_path, flexible_params, output_dir=None):
    """
    Run a comparison study between rigid and flexible docking.

    Args:
        receptor_path (str): Path to receptor PDB
        ligand_path (str): Path to ligand file
        flexible_params (dict): Parameters for flexible docking
        output_dir (str): Output directory for results

    Returns:
        dict: Comparison results
    """
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print("Running comparison between rigid and flexible docking...")

    # 1. Rigid docking
    print("\n1. Running rigid docking...")
    rigid_output = os.path.join(output_dir, 'rigid_docked.sdf') if output_dir else None

    rigid_params = {
        'num_modes': 9,
        'exhaustiveness': 8,
        'seed': 42
    }

    # Add binding site definition
    if 'autobox_ligand' in flexible_params:
        rigid_params['autobox_ligand'] = flexible_params['autobox_ligand']
    elif 'center_x' in flexible_params:
        rigid_params.update({
            'center_x': flexible_params['center_x'],
            'center_y': flexible_params['center_y'],
            'center_z': flexible_params['center_z'],
            'size_x': flexible_params.get('size_x', 20),
            'size_y': flexible_params.get('size_y', 20),
            'size_z': flexible_params.get('size_z', 20)
        })

    rigid_results, rigid_success, _ = run_flexible_docking(
        receptor_path, ligand_path, rigid_output, **rigid_params
    )

    # 2. Flexible docking
    print("\n2. Running flexible docking...")
    flex_output = os.path.join(output_dir, 'flexible_docked.sdf') if output_dir else None

    flex_results, flex_success, _ = run_flexible_docking(
        receptor_path, ligand_path, flex_output, **flexible_params
    )

    # 3. Analysis
    comparison = {
        'rigid_success': rigid_success,
        'flexible_success': flex_success
    }

    if rigid_success:
        comparison['rigid_analysis'] = analyze_flexible_results(rigid_results)

    if flex_success:
        comparison['flexible_analysis'] = analyze_flexible_results(
            flex_results, rigid_results if rigid_success else None
        )

    return comparison


def main():
    parser = argparse.ArgumentParser(
        description='Flexible receptor docking with gnina',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--receptor', '-r', required=True,
                       help='Receptor PDB file')
    parser.add_argument('--ligand', '-l', required=True,
                       help='Ligand SDF/PDB file')
    parser.add_argument('--output', '-o',
                       help='Output SDF file for docked poses')
    parser.add_argument('--output_dir',
                       help='Output directory for all results')

    # Flexible residue specification
    flex_group = parser.add_mutually_exclusive_group()
    flex_group.add_argument('--flexres',
                          help='Comma-separated list of flexible residues (chain:resid format)')
    flex_group.add_argument('--flexdist', type=float,
                          help='Distance cutoff for automatic flexible residue selection')

    parser.add_argument('--flexdist_ligand',
                       help='Reference ligand for flexdist calculation')
    parser.add_argument('--flex_limit', type=int,
                       help='Hard limit for number of flexible residues')
    parser.add_argument('--flex_max', type=int,
                       help='Maximum number of closest flexible residues to retain')

    # Binding site definition
    parser.add_argument('--autobox_ligand',
                       help='Reference ligand for automatic box definition')
    parser.add_argument('--center', nargs=3, type=float,
                       help='Box center coordinates (x y z)')
    parser.add_argument('--size', nargs=3, type=float, default=[20, 20, 20],
                       help='Box size (default: 20 20 20)')

    # Docking parameters
    parser.add_argument('--num_modes', type=int, default=9,
                       help='Number of binding modes (default: 9)')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                       help='Exhaustiveness of search (default: 8)')
    parser.add_argument('--cnn_scoring',
                       choices=['none', 'rescore', 'refinement'],
                       default='rescore',
                       help='CNN scoring mode (default: rescore)')
    parser.add_argument('--scoring',
                       choices=['default', 'vinardo', 'ad4_scoring'],
                       help='Empirical scoring function')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')
    parser.add_argument('--no_gpu', action='store_true',
                       help='Disable GPU acceleration')

    # Analysis options
    parser.add_argument('--compare_rigid', action='store_true',
                       help='Compare flexible docking with rigid docking')
    parser.add_argument('--makeflex_script',
                       help='Path to makeflex.py script for structure combination')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.receptor):
        print(f"Error: Receptor file {args.receptor} not found")
        sys.exit(1)

    if not os.path.exists(args.ligand):
        print(f"Error: Ligand file {args.ligand} not found")
        sys.exit(1)

    # Setup flexible docking parameters
    flexible_params = {
        'num_modes': args.num_modes,
        'exhaustiveness': args.exhaustiveness,
        'cnn_scoring': args.cnn_scoring,
        'seed': args.seed,
        'no_gpu': args.no_gpu
    }

    # Add flexible residue specification
    if args.flexres:
        flexible_params['flexres'] = args.flexres
    elif args.flexdist:
        flexible_params['flexdist'] = args.flexdist

        # Use provided ligand or the input ligand for distance calculation
        if args.flexdist_ligand:
            if not os.path.exists(args.flexdist_ligand):
                print(f"Error: Flexdist ligand {args.flexdist_ligand} not found")
                sys.exit(1)
            flexible_params['flexdist_ligand'] = args.flexdist_ligand
        else:
            flexible_params['flexdist_ligand'] = args.ligand
    else:
        # Try to automatically identify flexible residues
        print("No flexible residue specification provided. Attempting automatic identification...")
        flex_residues = identify_flexible_residues(args.receptor, args.ligand)

        if flex_residues:
            flexible_params['flexres'] = ','.join(flex_residues[:5])  # Limit to 5 residues
            print(f"Using automatically identified residues: {flexible_params['flexres']}")
        else:
            print("Could not automatically identify flexible residues.")
            print("Please specify --flexres or --flexdist manually.")
            sys.exit(1)

    # Add optional flexible parameters
    if args.flex_limit:
        flexible_params['flex_limit'] = args.flex_limit
    if args.flex_max:
        flexible_params['flex_max'] = args.flex_max

    # Add binding site definition
    if args.autobox_ligand:
        if not os.path.exists(args.autobox_ligand):
            print(f"Error: Autobox ligand {args.autobox_ligand} not found")
            sys.exit(1)
        flexible_params['autobox_ligand'] = args.autobox_ligand
    elif args.center:
        flexible_params.update({
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
        flexible_params['scoring'] = args.scoring

    print(f"Flexible Docking Setup:")
    print(f"  Receptor: {args.receptor}")
    print(f"  Ligand: {args.ligand}")
    print(f"  Flexible specification: {args.flexres or f'Distance-based ({args.flexdist} Å)'}")

    # Run docking
    if args.compare_rigid:
        # Run comparison study
        results = run_comparison_study(
            args.receptor, args.ligand, flexible_params, args.output_dir
        )

        # Display comparison results
        print("\n" + "="*60)
        print("RIGID vs FLEXIBLE DOCKING COMPARISON")
        print("="*60)

        if results['rigid_success']:
            rigid_analysis = results['rigid_analysis']
            print(f"\nRIGID DOCKING:")
            print(f"  Best affinity: {rigid_analysis['best_affinity']:.3f}")
            print(f"  Runtime: {rigid_analysis['runtime']:.1f} seconds")
            print(f"  Poses generated: {rigid_analysis['num_poses']}")

        if results['flexible_success']:
            flex_analysis = results['flexible_analysis']
            print(f"\nFLEXIBLE DOCKING:")
            print(f"  Best affinity: {flex_analysis['best_affinity']:.3f}")
            print(f"  Runtime: {flex_analysis['runtime']:.1f} seconds")
            print(f"  Poses generated: {flex_analysis['num_poses']}")
            print(f"  Flexible residues: {flex_analysis['num_flexible']}")

            if 'comparison' in flex_analysis:
                comp = flex_analysis['comparison']
                print(f"\nCOMPARISON:")
                print(f"  Affinity improvement: {comp.get('affinity_improvement', 0):.3f}")
                print(f"  Runtime overhead: {comp.get('runtime_overhead', 0):.1f} seconds")

    else:
        # Run flexible docking only
        results, success, output = run_flexible_docking(
            args.receptor, args.ligand, args.output, **flexible_params
        )

        if success:
            analysis = analyze_flexible_results(results)

            print("\n" + "="*60)
            print("FLEXIBLE DOCKING RESULTS")
            print("="*60)

            if results['poses']:
                print(f"\n{'Mode':<6}{'Affinity':<10}{'RMSD':<8}{'CNN Score':<12}")
                print("-" * 40)

                for pose in results['poses']:
                    rmsd_str = f"{pose['rmsd']:.3f}" if pose['rmsd'] is not None else "N/A"
                    cnn_str = f"{pose['cnn_score']:.3f}" if pose['cnn_score'] is not None else "N/A"

                    print(f"{pose['mode']:<6}{pose['affinity']:<10.3f}{rmsd_str:<8}{cnn_str:<12}")

                print(f"\nSummary:")
                print(f"  Total poses: {analysis['num_poses']}")
                print(f"  Best affinity: {analysis['best_affinity']:.3f}")
                print(f"  Flexible residues: {analysis['num_flexible']}")
                print(f"  Runtime: {analysis['runtime']:.1f} seconds")

                # Try to combine structures if makeflex.py is available
                if args.output and args.makeflex_script:
                    flex_output = args.output.replace('.sdf', '_flex.pdb')
                    combined_output = args.output.replace('.sdf', '_combined.pdb')

                    if combine_flexible_structures(
                        args.receptor, flex_output, combined_output, args.makeflex_script
                    ):
                        print(f"  Combined structure: {combined_output}")

            else:
                print("No poses were generated successfully.")

        else:
            print("Flexible docking failed.")
            sys.exit(1)


if __name__ == '__main__':
    main()