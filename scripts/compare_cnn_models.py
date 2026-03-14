#!/usr/bin/env python3
"""
Script: compare_cnn_models.py
Description: Compare CNN model performance for protein-ligand scoring

Original Use Case: examples/use_case_4_cnn_model_comparison.py
Dependencies Removed: matplotlib/seaborn (optional plotting)

Usage:
    python scripts/compare_cnn_models.py --receptor <receptor.pdb> --ligand <ligand.sdf> --output <comparison.csv>

Example:
    python scripts/compare_cnn_models.py --receptor examples/data/184l_rec.pdb --ligand examples/data/184l_lig.sdf --output results/cnn_comparison.csv
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
import time
from pathlib import Path
from typing import Union, Optional, Dict, Any, List, Tuple

# Essential scientific packages
import pandas as pd
import numpy as np

# ==============================================================================
# Configuration
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
    "comparison": {
        "models": ["default", "fast", "dense"],
        "modes": ["score"],
        "iterations": 3,
        "include_timing": True
    },
    "benchmarking": {
        "warmup_runs": 1,
        "statistical_analysis": True,
        "correlation_analysis": False  # Simplified for MCP
    },
    "output": {
        "format": "csv",
        "detailed_results": True,
        "summary_table": True
    },
    "gnina": {
        "executable": _resolve_gnina_executable(),
        "timeout": 900,
        "verbose": False
    }
}

# ==============================================================================
# Inlined Utility Functions
# ==============================================================================
def parse_cnn_scores(output_text: str) -> Dict[str, float]:
    """Parse CNN scores from gnina output."""
    output = output_text.decode() if isinstance(output_text, bytes) else output_text
    scores = {}

    # Extract different score types
    patterns = {
        'affinity': r'Affinity:\s+(\S+)',
        'cnn_score': r'CNNscore:\s+(\S+)',
        'cnn_affinity': r'CNNaffinity:\s+(\S+)'
    }

    for score_type, pattern in patterns.items():
        match = re.search(pattern, output)
        if match:
            try:
                scores[score_type] = float(match.group(1))
            except ValueError:
                scores[score_type] = None
        else:
            scores[score_type] = None

    return scores


def benchmark_cnn_model(receptor_path: str, ligand_path: str, cnn_model: str,
                       mode: str, iterations: int, config: Dict[str, Any]) -> Dict[str, Any]:
    """Benchmark a specific CNN model."""
    results = {
        'model': cnn_model,
        'mode': mode,
        'iterations': iterations,
        'successful_runs': 0,
        'failed_runs': 0,
        'scores': [],
        'runtimes': [],
        'errors': []
    }

    for i in range(iterations):
        try:
            start_time = time.time()

            # Build command
            cmd = [
                config['gnina']['executable'],
                '-r', str(receptor_path),
                '-l', str(ligand_path)
            ]

            if mode == 'score':
                cmd.append('--score_only')
            elif mode == 'dock':
                cmd.extend(['--num_modes', '3', '--exhaustiveness', '4'])

            # Add CNN model
            if cnn_model != 'default':
                cmd.extend(['--cnn', cnn_model])

            # Add CNN scoring
            cmd.extend(['--cnn_scoring', 'all'])

            # Run gnina
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True,
                timeout=config['gnina']['timeout'],
                env=_gnina_subprocess_env()
            )

            runtime = time.time() - start_time
            output = result.stdout + result.stderr

            # Parse scores
            scores = parse_cnn_scores(output)
            scores['runtime'] = runtime
            scores['iteration'] = i + 1

            results['scores'].append(scores)
            results['runtimes'].append(runtime)
            results['successful_runs'] += 1

        except Exception as e:
            results['failed_runs'] += 1
            results['errors'].append(str(e))

    return results


def analyze_benchmark_results(results_list: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze and compare benchmark results across models."""
    analysis = {
        'model_comparison': [],
        'runtime_stats': {},
        'score_stats': {},
        'recommendations': []
    }

    for result in results_list:
        if result['successful_runs'] == 0:
            continue

        model_name = result['model']
        scores = result['scores']
        runtimes = result['runtimes']

        # Runtime statistics
        analysis['runtime_stats'][model_name] = {
            'mean_runtime': np.mean(runtimes),
            'std_runtime': np.std(runtimes),
            'min_runtime': np.min(runtimes),
            'max_runtime': np.max(runtimes)
        }

        # Score statistics
        affinities = [s['affinity'] for s in scores if s['affinity'] is not None]
        cnn_scores = [s['cnn_score'] for s in scores if s['cnn_score'] is not None]
        cnn_affinities = [s['cnn_affinity'] for s in scores if s['cnn_affinity'] is not None]

        analysis['score_stats'][model_name] = {
            'mean_affinity': np.mean(affinities) if affinities else None,
            'std_affinity': np.std(affinities) if affinities else None,
            'mean_cnn_score': np.mean(cnn_scores) if cnn_scores else None,
            'std_cnn_score': np.std(cnn_scores) if cnn_scores else None,
            'mean_cnn_affinity': np.mean(cnn_affinities) if cnn_affinities else None
        }

        # Model comparison summary
        analysis['model_comparison'].append({
            'model': model_name,
            'success_rate': result['successful_runs'] / (result['successful_runs'] + result['failed_runs']),
            'mean_runtime': np.mean(runtimes),
            'mean_affinity': np.mean(affinities) if affinities else None,
            'reliability': result['successful_runs'] / len(runtimes) if runtimes else 0
        })

    # Generate recommendations
    if analysis['model_comparison']:
        # Fastest model
        fastest = min(analysis['model_comparison'], key=lambda x: x['mean_runtime'])
        analysis['recommendations'].append(f"Fastest model: {fastest['model']} ({fastest['mean_runtime']:.2f}s)")

        # Most reliable model
        most_reliable = max(analysis['model_comparison'], key=lambda x: x['reliability'])
        analysis['recommendations'].append(f"Most reliable: {most_reliable['model']} ({most_reliable['reliability']*100:.1f}% success)")

        # Best affinity predictor (if available)
        valid_affinities = [m for m in analysis['model_comparison'] if m['mean_affinity'] is not None]
        if valid_affinities:
            best_affinity = min(valid_affinities, key=lambda x: x['mean_affinity'])
            analysis['recommendations'].append(f"Best affinity prediction: {best_affinity['model']} ({best_affinity['mean_affinity']:.3f} kcal/mol)")

    return analysis


# ==============================================================================
# Core Function
# ==============================================================================
def run_cnn_comparison(
    receptor_file: Union[str, Path],
    ligand_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    models: Optional[List[str]] = None,
    modes: Optional[List[str]] = None,
    iterations: int = 3,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Compare CNN model performance for protein-ligand scoring.

    Args:
        receptor_file: Path to receptor PDB file
        ligand_file: Path to ligand file
        output_file: Path to save comparison results
        models: List of CNN models to compare
        modes: List of modes to test ('score', 'dock')
        iterations: Number of iterations per model
        config: Configuration dict
        **kwargs: Override specific config parameters

    Returns:
        Dict containing comparison results and analysis

    Example:
        >>> result = run_cnn_comparison("receptor.pdb", "ligand.sdf",
        ...                           models=["default", "fast"])
        >>> print(result['analysis']['recommendations'])
    """
    # Setup
    receptor_file = Path(receptor_file)
    ligand_file = Path(ligand_file)
    config = _deep_merge(DEFAULT_CONFIG, config or {})
    if kwargs:
        config = _deep_merge(config, kwargs)

    # Override parameters
    if models:
        config['comparison']['models'] = models
    if modes:
        config['comparison']['modes'] = modes
    config['comparison']['iterations'] = iterations

    # Validate inputs
    if not receptor_file.exists():
        raise FileNotFoundError(f"Receptor file not found: {receptor_file}")
    if not ligand_file.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")

    print(f"Comparing {len(config['comparison']['models'])} CNN models...")

    # Run benchmarks
    all_results = []
    total_tests = len(config['comparison']['models']) * len(config['comparison']['modes'])
    current_test = 0

    for model in config['comparison']['models']:
        for mode in config['comparison']['modes']:
            current_test += 1
            print(f"Testing {model} in {mode} mode ({current_test}/{total_tests})...")

            # Run warmup if configured
            if config['benchmarking']['warmup_runs'] > 0:
                for _ in range(config['benchmarking']['warmup_runs']):
                    try:
                        benchmark_cnn_model(receptor_file, ligand_file, model, mode, 1, config)
                    except:
                        pass  # Ignore warmup failures

            # Run actual benchmark
            result = benchmark_cnn_model(
                receptor_file, ligand_file, model, mode,
                config['comparison']['iterations'], config
            )
            all_results.append(result)

    # Analyze results
    analysis = analyze_benchmark_results(all_results)

    # Prepare output data
    detailed_results = []
    for result in all_results:
        for score_data in result['scores']:
            row = {
                'model': result['model'],
                'mode': result['mode'],
                'iteration': score_data['iteration'],
                'runtime': score_data['runtime'],
                'affinity': score_data['affinity'],
                'cnn_score': score_data['cnn_score'],
                'cnn_affinity': score_data['cnn_affinity']
            }
            detailed_results.append(row)

    # Save results
    output_path = None
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save detailed results
        df = pd.DataFrame(detailed_results)
        df.to_csv(output_path, index=False, float_format='%.4f')

        # Save summary
        if config['output']['summary_table']:
            summary_path = output_path.parent / f"{output_path.stem}_summary.csv"
            summary_df = pd.DataFrame(analysis['model_comparison'])
            summary_df.to_csv(summary_path, index=False, float_format='%.4f')

    return {
        "detailed_results": detailed_results,
        "benchmark_results": all_results,
        "analysis": analysis,
        "output_file": str(output_path) if output_path else None,
        "metadata": {
            "receptor_file": str(receptor_file),
            "ligand_file": str(ligand_file),
            "total_tests": total_tests,
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
    parser.add_argument('--ligand', '-l', required=True, help='Ligand file')
    parser.add_argument('--output', '-o', help='Output CSV file')
    parser.add_argument('--models', nargs='+',
                       default=['default', 'fast'],
                       choices=['default', 'fast', 'dense', 'default1.0', 'general_default2018'],
                       help='CNN models to compare')
    parser.add_argument('--modes', nargs='+',
                       default=['score'],
                       choices=['score', 'dock'],
                       help='Test modes')
    parser.add_argument('--iterations', type=int, default=3,
                       help='Iterations per model')
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
        "gnina": {
            "verbose": args.verbose
        }
    }

    if config:
        config['gnina'].update(config_overrides['gnina'])
    else:
        config = config_overrides

    try:
        output_path = Path(args.output) if args.output else None

        # If output is .json, use a sibling .csv for results and write JSON
        if output_path and output_path.suffix.lower() == '.json':
            csv_output = str(output_path.with_suffix('.csv'))
        else:
            csv_output = args.output

        # Run CNN comparison
        result = run_cnn_comparison(
            receptor_file=args.receptor,
            ligand_file=args.ligand,
            output_file=csv_output,
            models=args.models,
            modes=args.modes,
            iterations=args.iterations,
            config=config
        )

        # Print results
        print("="*60)
        print("CNN MODEL COMPARISON RESULTS")
        print("="*60)

        # Summary table
        if result['analysis']['model_comparison']:
            df = pd.DataFrame(result['analysis']['model_comparison'])
            print(df.to_string(index=False, float_format='%.3f'))

        print("\n" + "="*60)
        print("RECOMMENDATIONS")
        print("="*60)
        for rec in result['analysis']['recommendations']:
            print(f"• {rec}")

        if result['output_file']:
            print(f"\nDetailed results saved to: {result['output_file']}")

        # Write JSON results if requested
        if output_path and output_path.suffix.lower() == '.json':
            result.pop('metadata', None)
            with open(output_path, 'w') as f:
                json.dump(result, f, indent=2)

        return result

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()