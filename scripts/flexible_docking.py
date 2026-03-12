#!/usr/bin/env python3
"""
Script: flexible_docking.py
Description: Flexible receptor docking with gnina

Original Use Case: examples/use_case_5_flexible_docking.py
Dependencies Removed: makeflex.py integration (simplified)

Usage:
    python scripts/flexible_docking.py --receptor <receptor.pdb> --ligand <ligand.sdf> --flexdist 3.5 --output <output.sdf>

Example:
    python scripts/flexible_docking.py --receptor examples/data/184l_rec.pdb --ligand examples/data/184l_lig.sdf --flexdist 3.5 --output results/flex_docked.sdf
"""

import argparse
import subprocess
import sys
import os
import json
from pathlib import Path
from typing import Union, Optional, Dict, Any, List

import pandas as pd

DEFAULT_CONFIG = {
    "flexibility": {
        "flexdist": 3.5,
        "flexres": None,
        "flex_residues": [],
        "flexdist_ligand": None
    },
    "docking": {
        "num_modes": 9,
        "exhaustiveness": 8,
        "autobox_add": 4,
        "compare_rigid": False
    },
    "output": {
        "format": "sdf",
        "flexible_receptor_pdb": False,
        "comparison_report": True
    },
    "gnina": {
        "executable": "gnina",
        "timeout": 2400,
        "gpu": True,
        "verbose": False
    }
}

def run_flexible_docking(
    receptor_file: Union[str, Path],
    ligand_file: Union[str, Path],
    output_file: Union[str, Path],
    flexdist: float = 3.5,
    flexdist_ligand: Optional[Union[str, Path]] = None,
    flexres: Optional[str] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform flexible receptor docking using gnina.

    Args:
        receptor_file: Path to receptor PDB file
        ligand_file: Path to ligand file
        output_file: Path to save docked poses
        flexdist: Distance cutoff for flexible residues
        flexdist_ligand: Reference ligand for flexdist calculation
        flexres: Specific residues to make flexible (e.g., "A:120,A:145")
        config: Configuration dict
        **kwargs: Override specific config parameters

    Returns:
        Dict containing docking results and analysis
    """
    # Setup
    receptor_file = Path(receptor_file)
    ligand_file = Path(ligand_file)
    output_file = Path(output_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Override parameters
    if flexdist:
        config['flexibility']['flexdist'] = flexdist
    if flexres:
        config['flexibility']['flexres'] = flexres
    if flexdist_ligand:
        config['flexibility']['flexdist_ligand'] = flexdist_ligand

    # Validate inputs
    if not receptor_file.exists():
        raise FileNotFoundError(f"Receptor file not found: {receptor_file}")
    if not ligand_file.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")

    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Build flexible docking command
    cmd = [
        config['gnina']['executable'],
        '-r', str(receptor_file),
        '-l', str(ligand_file),
        '-o', str(output_file),
        '--num_modes', str(config['docking']['num_modes']),
        '--exhaustiveness', str(config['docking']['exhaustiveness'])
    ]

    # Add flexibility options
    if config['flexibility']['flexres']:
        cmd.extend(['--flexres', config['flexibility']['flexres']])
    elif config['flexibility']['flexdist_ligand']:
        cmd.extend(['--flexdist_ligand', str(config['flexibility']['flexdist_ligand'])])
        cmd.extend(['--flexdist', str(config['flexibility']['flexdist'])])
    else:
        # Use flexdist with the docking ligand
        cmd.extend(['--flexdist_ligand', str(ligand_file)])
        cmd.extend(['--flexdist', str(config['flexibility']['flexdist'])])

    # Add autobox
    cmd.extend(['--autobox_ligand', str(ligand_file)])
    cmd.extend(['--autobox_add', str(config['docking']['autobox_add'])])

    if config['gnina']['gpu']:
        cmd.append('--gpu')

    try:
        # Run flexible docking
        if config['gnina']['verbose']:
            print(f"Running flexible docking: {' '.join(cmd)}")

        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True,
            timeout=config['gnina']['timeout']
        )
        output = result.stdout + result.stderr

        # Parse results (simplified)
        poses = []
        lines = output.split('\n')
        for line in lines:
            if line.strip().startswith(('1', '2', '3', '4', '5', '6', '7', '8', '9')):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        pose_num = int(parts[0])
                        affinity = float(parts[1])
                        poses.append({"pose": pose_num, "affinity": affinity})
                    except (ValueError, IndexError):
                        continue

        # Compare with rigid docking if requested
        rigid_poses = []
        if config['docking']['compare_rigid']:
            rigid_cmd = [c for c in cmd if not c.startswith('--flex')]
            rigid_output_file = output_file.parent / f"{output_file.stem}_rigid.sdf"
            rigid_cmd[-3] = str(rigid_output_file)  # Replace output file

            try:
                rigid_result = subprocess.run(
                    rigid_cmd, capture_output=True, text=True, check=True,
                    timeout=config['gnina']['timeout']
                )
                # Parse rigid results (simplified)
                rigid_lines = rigid_result.stdout.split('\n') + rigid_result.stderr.split('\n')
                for line in rigid_lines:
                    if line.strip().startswith(('1', '2', '3', '4', '5', '6', '7', '8', '9')):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                pose_num = int(parts[0])
                                affinity = float(parts[1])
                                rigid_poses.append({"pose": pose_num, "affinity": affinity, "type": "rigid"})
                            except (ValueError, IndexError):
                                continue
            except Exception as e:
                print(f"Warning: Rigid comparison failed: {e}")

        # Generate analysis
        analysis = {
            "flexible_poses": len(poses),
            "best_flexible_affinity": min([p['affinity'] for p in poses]) if poses else None,
            "rigid_poses": len(rigid_poses),
            "best_rigid_affinity": min([p['affinity'] for p in rigid_poses]) if rigid_poses else None
        }

        if analysis['best_flexible_affinity'] and analysis['best_rigid_affinity']:
            analysis['improvement'] = analysis['best_rigid_affinity'] - analysis['best_flexible_affinity']

        return {
            "flexible_poses": poses,
            "rigid_poses": rigid_poses,
            "analysis": analysis,
            "output_file": str(output_file),
            "metadata": {
                "receptor_file": str(receptor_file),
                "ligand_file": str(ligand_file),
                "flexdist": config['flexibility']['flexdist'],
                "config": config
            }
        }

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Flexible docking timed out after {config['gnina']['timeout']} seconds")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Flexible docking failed: {e.stderr}")

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--receptor', '-r', required=True, help='Receptor PDB file')
    parser.add_argument('--ligand', '-l', required=True, help='Ligand file')
    parser.add_argument('--output', '-o', required=True, help='Output SDF file')
    parser.add_argument('--flexdist', type=float, default=3.5, help='Flexibility distance cutoff')
    parser.add_argument('--flexdist_ligand', help='Reference ligand for flexdist')
    parser.add_argument('--flexres', help='Specific residues to make flexible')
    parser.add_argument('--compare_rigid', action='store_true', help='Compare with rigid docking')
    parser.add_argument('--config', '-c', help='Config file (JSON)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    try:
        result = run_flexible_docking(
            receptor_file=args.receptor,
            ligand_file=args.ligand,
            output_file=args.output,
            flexdist=args.flexdist,
            flexdist_ligand=args.flexdist_ligand,
            flexres=args.flexres,
            config=config,
            docking={'compare_rigid': args.compare_rigid},
            gnina={'verbose': args.verbose}
        )

        print("="*60)
        print("FLEXIBLE DOCKING RESULTS")
        print("="*60)

        if result['flexible_poses']:
            print(f"Flexible poses: {len(result['flexible_poses'])}")
            print(f"Best flexible affinity: {result['analysis']['best_flexible_affinity']:.3f} kcal/mol")

        if result['rigid_poses']:
            print(f"Rigid poses: {len(result['rigid_poses'])}")
            print(f"Best rigid affinity: {result['analysis']['best_rigid_affinity']:.3f} kcal/mol")
            if result['analysis'].get('improvement'):
                print(f"Improvement: {result['analysis']['improvement']:.3f} kcal/mol")

        print(f"\nOutput saved to: {result['output_file']}")

        return result

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()