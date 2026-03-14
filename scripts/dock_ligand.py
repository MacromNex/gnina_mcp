#!/usr/bin/env python3
"""
Script: dock_ligand.py
Description: Perform molecular docking with gnina CNN-enhanced scoring

Original Use Case: examples/use_case_2_standard_docking.py
Dependencies Removed: None (already minimal)

Usage:
    python scripts/dock_ligand.py --receptor <receptor.pdb> --ligand <ligand.sdf> --output <output.sdf>

Example:
    python scripts/dock_ligand.py --receptor examples/data/184l_rec.pdb --ligand examples/data/184l_lig.sdf --autobox_ligand examples/data/184l_lig.sdf --output results/docked.sdf
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import shutil
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
def _resolve_gnina_executable():
    """Find gnina executable, falling back to env/bin or mock_gnina.py."""
    if shutil.which('gnina'):
        return 'gnina'
    # Check the conda env bundled with the project
    env_path = Path(__file__).resolve().parent.parent / 'env' / 'bin' / 'gnina'
    if env_path.exists():
        return str(env_path)
    # Fall back to mock_gnina.py for testing
    mock_path = Path(__file__).resolve().parent.parent / 'mock_gnina.py'
    if mock_path.exists():
        return str(mock_path)
    return 'gnina'

def _gnina_subprocess_env():
    """Get environment dict with LD_LIBRARY_PATH set for gnina."""
    env = os.environ.copy()
    env_lib = Path(__file__).resolve().parent.parent / 'env' / 'lib'
    if env_lib.exists():
        ld_path = env.get('LD_LIBRARY_PATH', '')
        if str(env_lib) not in ld_path:
            env['LD_LIBRARY_PATH'] = f"{env_lib}:{ld_path}" if ld_path else str(env_lib)
    return env

def _deep_merge(base: dict, override: dict) -> dict:
    """Deep merge override into base, returning a new dict."""
    result = base.copy()
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result

DEFAULT_CONFIG = {
    "docking": {
        "num_modes": 9,
        "exhaustiveness": 8,
        "seed": 0,
        "autobox_add": 4,
        "spacing": 0.375
    },
    "scoring": {
        "cnn_scoring": "rescore",
        "cnn_model": "default",
        "empirical_only": False
    },
    "binding_site": {
        "autobox_ligand": None,
        "center": None,
        "size": [20, 20, 20]
    },
    "output": {
        "format": "sdf",
        "include_hydrogens": True,
        "include_scores": True
    },
    "gnina": {
        "executable": _resolve_gnina_executable(),
        "timeout": 1800,
        "verbose": False
    }
}

# ==============================================================================
# Inlined Utility Functions
# ==============================================================================
def parse_docking_output(output_text: str) -> List[Dict[str, Any]]:
    """Parse gnina docking output to extract pose information.

    Gnina output format (docking mode with CNN rescore):
        mode |  affinity  |  intramol  |    CNN     |   CNN
             | (kcal/mol) | (kcal/mol) | pose score | affinity
        -----+------------+------------+------------+----------
            1       -6.87       -0.31       0.9342      4.874

    Gnina output format (docking mode without CNN):
        mode |   affinity | dist from best mode
             | (kcal/mol) | rmsd l.b.| rmsd u.b.
        -----+------------+----------+----------
            1       -6.87      0.000      0.000
    """
    output = output_text.decode() if isinstance(output_text, bytes) else output_text
    poses = []

    # Match pose lines: leading whitespace + mode number + numeric columns
    # This avoids matching progress bar "0%   10   20   30..." lines
    pose_pattern = r'^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)(?:\s+(\d+\.\d+))?'
    for line in output.splitlines():
        match = re.match(pose_pattern, line)
        if match:
            try:
                pose_num = int(match.group(1))
                affinity = float(match.group(2))
                col3 = float(match.group(3))
                col4 = float(match.group(4))
                col5 = float(match.group(5)) if match.group(5) else None

                pose = {
                    "pose": pose_num,
                    "affinity": affinity,
                }

                if col5 is not None:
                    # 5-column CNN format: mode, affinity, intramol, cnn_score, cnn_affinity
                    pose["intramol"] = col3
                    pose["cnn_score"] = col4
                    pose["cnn_affinity"] = col5
                else:
                    # 4-column Vina format: mode, affinity, rmsd_lb, rmsd_ub
                    pose["rmsd_lb"] = col3
                    pose["rmsd_ub"] = col4

                poses.append(pose)
            except (ValueError, IndexError):
                continue

    return poses


def validate_docking_inputs(receptor_path: Union[str, Path], ligand_path: Union[str, Path],
                          autobox_ligand: Optional[Union[str, Path]] = None,
                          center: Optional[List[float]] = None) -> bool:
    """Validate docking input parameters."""
    receptor_path = Path(receptor_path)
    ligand_path = Path(ligand_path)

    if not receptor_path.exists():
        raise FileNotFoundError(f"Receptor file not found: {receptor_path}")

    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_path}")

    # Validate binding site specification
    if autobox_ligand is not None:
        autobox_path = Path(autobox_ligand)
        if not autobox_path.exists():
            raise FileNotFoundError(f"Autobox ligand file not found: {autobox_path}")
    elif center is None:
        raise ValueError("Either autobox_ligand or center coordinates must be specified")

    if center is not None and len(center) != 3:
        raise ValueError("Center coordinates must be [x, y, z]")

    return True


def run_docking_command(receptor_path: str, ligand_path: str, output_path: str,
                       config: Dict[str, Any]) -> Tuple[List[Dict], str]:
    """Execute gnina docking command with error handling."""
    # Build gnina command
    cmd = [
        config['gnina']['executable'],
        '-r', str(receptor_path),
        '-l', str(ligand_path),
        '-o', str(output_path),
        '--num_modes', str(config['docking']['num_modes']),
        '--exhaustiveness', str(config['docking']['exhaustiveness'])
    ]

    # Add seed if specified
    if config['docking'].get('seed') is not None:
        cmd.extend(['--seed', str(config['docking']['seed'])])

    # Add binding site specification
    if config['binding_site']['autobox_ligand']:
        cmd.extend(['--autobox_ligand', str(config['binding_site']['autobox_ligand'])])
        cmd.extend(['--autobox_add', str(config['docking']['autobox_add'])])
    elif config['binding_site']['center']:
        center = config['binding_site']['center']
        size = config['binding_site']['size']
        cmd.extend(['--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2])])
        cmd.extend(['--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2])])

    # Add CNN scoring options
    if config['scoring']['cnn_scoring'] != 'none':
        cmd.extend(['--cnn_scoring', config['scoring']['cnn_scoring']])
    if config['scoring']['cnn_model'] != 'default':
        cmd.extend(['--cnn', config['scoring']['cnn_model']])

    try:
        # Run gnina docking
        if config['gnina']['verbose']:
            print(f"Running command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True,
            timeout=config['gnina']['timeout'],
            env=_gnina_subprocess_env()
        )
        output = result.stdout + result.stderr

        # Parse poses
        poses = parse_docking_output(output)

        return poses, output

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Docking command timed out after {config['gnina']['timeout']} seconds")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Docking command failed: {e.stderr}")
    except FileNotFoundError:
        raise RuntimeError(f"Gnina executable not found: {config['gnina']['executable']}")


def analyze_poses(poses: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze docking poses and generate statistics."""
    if not poses:
        return {"error": "No poses generated"}

    affinities = [pose['affinity'] for pose in poses if pose['affinity'] is not None]
    cnn_affinities = [pose.get('cnn_affinity') for pose in poses if pose.get('cnn_affinity') is not None]

    analysis = {
        "total_poses": len(poses),
        "best_affinity": min(affinities) if affinities else None,
        "worst_affinity": max(affinities) if affinities else None,
        "mean_affinity": sum(affinities) / len(affinities) if affinities else None,
        "affinity_range": max(affinities) - min(affinities) if len(affinities) > 1 else 0
    }

    if cnn_affinities:
        analysis.update({
            "best_cnn_affinity": max(cnn_affinities),
            "worst_cnn_affinity": min(cnn_affinities),
            "mean_cnn_affinity": sum(cnn_affinities) / len(cnn_affinities)
        })

    return analysis


# ==============================================================================
# Core Function
# ==============================================================================
def run_molecular_docking(
    receptor_file: Union[str, Path],
    ligand_file: Union[str, Path],
    output_file: Union[str, Path],
    autobox_ligand: Optional[Union[str, Path]] = None,
    center: Optional[List[float]] = None,
    size: Optional[List[float]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform molecular docking using gnina with CNN-enhanced scoring.

    Args:
        receptor_file: Path to receptor PDB file
        ligand_file: Path to ligand file (SDF, PDB, etc.)
        output_file: Path to save docked poses (SDF format)
        autobox_ligand: Reference ligand for automatic binding site detection
        center: Binding site center coordinates [x, y, z]
        size: Binding site box size [x, y, z] (default from config)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - poses: List of pose information
            - analysis: Statistical analysis of poses
            - output_file: Path to output file
            - metadata: Execution metadata

    Example:
        >>> result = run_molecular_docking("receptor.pdb", "ligand.sdf", "docked.sdf",
        ...                               autobox_ligand="ref_ligand.sdf")
        >>> print(f"Generated {len(result['poses'])} poses")
    """
    # Setup
    receptor_file = Path(receptor_file)
    ligand_file = Path(ligand_file)
    output_file = Path(output_file)
    config = _deep_merge(DEFAULT_CONFIG, config or {})
    if kwargs:
        config = _deep_merge(config, kwargs)

    # Override binding site parameters if provided
    if autobox_ligand:
        config['binding_site']['autobox_ligand'] = autobox_ligand
    if center:
        config['binding_site']['center'] = center
    if size:
        config['binding_site']['size'] = size

    # Validate inputs
    validate_docking_inputs(receptor_file, ligand_file,
                          config['binding_site']['autobox_ligand'],
                          config['binding_site']['center'])

    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Run docking
    poses, raw_output = run_docking_command(receptor_file, ligand_file, output_file, config)

    # Analyze results
    analysis = analyze_poses(poses)

    return {
        "poses": poses,
        "analysis": analysis,
        "output_file": str(output_file),
        "metadata": {
            "receptor_file": str(receptor_file),
            "ligand_file": str(ligand_file),
            "config": config,
            "raw_output": raw_output if config['gnina']['verbose'] else None
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
    parser.add_argument('--output', '-o', required=True, help='Output SDF file for docked poses')
    parser.add_argument('--autobox_ligand', help='Reference ligand for autobox')
    parser.add_argument('--center', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Binding site center coordinates')
    parser.add_argument('--size', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Binding site box size')
    parser.add_argument('--num_modes', type=int, default=9,
                       help='Number of poses to generate')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                       help='Exhaustiveness of the global search')
    parser.add_argument('--cnn_scoring', default='rescore',
                       choices=['none', 'rescore', 'refinement', 'all'],
                       help='CNN scoring mode')
    parser.add_argument('--config', '-c', help='Config file (JSON)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    # Override config with command line arguments
    config_overrides = {
        "docking": {
            "num_modes": args.num_modes,
            "exhaustiveness": args.exhaustiveness
        },
        "scoring": {
            "cnn_scoring": args.cnn_scoring
        },
        "gnina": {
            "verbose": args.verbose
        }
    }

    # Merge configs
    if config:
        # Deep merge
        for section, params in config_overrides.items():
            if section not in config:
                config[section] = {}
            config[section].update(params)
    else:
        config = config_overrides

    try:
        output_path = Path(args.output)

        # If output is .json, use a sibling .sdf for gnina and write JSON results
        if output_path.suffix.lower() == '.json':
            sdf_output = str(output_path.with_suffix('.sdf'))
        else:
            sdf_output = args.output

        # Run docking
        result = run_molecular_docking(
            receptor_file=args.receptor,
            ligand_file=args.ligand,
            output_file=sdf_output,
            autobox_ligand=args.autobox_ligand,
            center=args.center,
            size=args.size,
            config=config
        )

        # Print results
        print("="*60)
        print("MOLECULAR DOCKING RESULTS")
        print("="*60)

        poses = result['poses']
        if poses:
            # Create DataFrame for nice display
            df = pd.DataFrame(poses)
            print(df.to_string(index=False, float_format='%.3f'))

            print(f"\nOutput file: {result['output_file']}")

            # Print analysis
            analysis = result['analysis']
            if 'error' not in analysis:
                print("\n" + "="*60)
                print("POSE ANALYSIS")
                print("="*60)
                print(f"Total poses: {analysis['total_poses']}")
                if analysis['best_affinity'] is not None:
                    print(f"Best affinity: {analysis['best_affinity']:.3f} kcal/mol")
                    print(f"Mean affinity: {analysis['mean_affinity']:.3f} kcal/mol")
                    print(f"Affinity range: {analysis['affinity_range']:.3f} kcal/mol")

            # Write JSON results if requested
            if output_path.suffix.lower() == '.json':
                # Remove non-serializable config from metadata
                result['metadata'].pop('raw_output', None)
                result['metadata'].pop('config', None)
                with open(output_path, 'w') as f:
                    json.dump(result, f, indent=2)
                print(f"JSON results saved to: {output_path}")
        else:
            print("No poses were generated successfully")
            sys.exit(1)

        return result

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()