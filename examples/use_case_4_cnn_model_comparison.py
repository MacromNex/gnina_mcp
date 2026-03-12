#!/usr/bin/env python3
"""
Use Case 4: CNN Model Comparison with Gnina
============================================

This script demonstrates how to compare different CNN models available in gnina
for protein-ligand scoring and binding affinity prediction. It evaluates multiple
CNN architectures on the same complexes to understand their performance characteristics.

Use case: Benchmark different CNN models for your specific system to choose the
best model for accuracy vs speed tradeoffs in cyclic peptide design.

Requirements:
- gnina executable installed and in PATH
- Receptor PDB file
- Ligand SDF file(s)
- Optional: Reference ligand for binding site definition

Example Usage:
    python use_case_4_cnn_model_comparison.py --receptor data/184l_rec.pdb --ligand data/184l_lig.sdf
    python use_case_4_cnn_model_comparison.py --receptor data/3rod_rec.pdb --ligand data/3rod_lig.pdb --models default fast dense --mode both
"""

import argparse
import subprocess
import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import time
import tempfile


def get_available_cnn_models():
    """
    Get list of available CNN models in gnina.

    Returns:
        dict: Dictionary mapping model names to descriptions
    """
    return {
        'default': 'Default ensemble (dense_1_3, dense_1_3_PT_KD_3, crossdock_default2018_KD_4)',
        'default1.0': 'GNINA 1.0 default ensemble (legacy)',
        'fast': 'Fast CNN model for rapid screening',
        'dense': 'Dense CNN architecture',
        'dense_1_3': 'Dense CNN variant 1.3',
        'dense_1_3_1': 'Dense CNN variant 1.3.1',
        'dense_1_3_2': 'Dense CNN variant 1.3.2',
        'dense_1_3_3': 'Dense CNN variant 1.3.3',
        'crossdock_default2018': 'CrossDocked 2018 default model',
        'crossdock_default2018_1_3': 'CrossDocked 2018 model 1.3',
        'general_default2018': 'General 2018 model',
        'redock_default2018': 'Redocking 2018 model',
        'all_default_to_default_1_3_1': 'All-dataset trained model 1.3.1'
    }


def parse_gnina_output_detailed(output_text):
    """
    Parse gnina output to extract detailed scoring information.

    Args:
        output_text (str): Raw output from gnina

    Returns:
        dict: Detailed scoring information
    """
    output = output_text.decode() if isinstance(output_text, bytes) else output_text

    # Initialize results
    results = {
        'affinity': None,
        'cnn_score': None,
        'cnn_affinity': None,
        'poses': [],
        'runtime': None
    }

    # Extract runtime information
    runtime_match = re.search(r'Refining results.*?(\d+\.?\d*)\s*seconds?', output, re.IGNORECASE)
    if runtime_match:
        results['runtime'] = float(runtime_match.group(1))

    # Extract single scores (for score-only mode)
    aff_match = re.search(r'Affinity:\s+(\S+)', output)
    if aff_match:
        results['affinity'] = float(aff_match.group(1))

    cnn_score_match = re.search(r'CNNscore:\s+(\S+)', output)
    if cnn_score_match:
        results['cnn_score'] = float(cnn_score_match.group(1))

    cnn_aff_match = re.search(r'CNNaffinity:\s+(\S+)', output)
    if cnn_aff_match:
        results['cnn_affinity'] = float(cnn_aff_match.group(1))

    # Extract pose table (for docking mode)
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

                    results['poses'].append({
                        'mode': mode,
                        'affinity': affinity,
                        'rmsd': rmsd,
                        'cnn_score': cnn_score,
                        'cnn_affinity': cnn_affinity
                    })

                    # Use first pose as representative
                    if mode == 1:
                        if not results['affinity']:
                            results['affinity'] = affinity
                        if not results['cnn_score']:
                            results['cnn_score'] = cnn_score
                        if not results['cnn_affinity']:
                            results['cnn_affinity'] = cnn_affinity

                except ValueError:
                    continue
            elif line.strip() == '':
                break

    return results


def run_cnn_comparison_single(receptor_path, ligand_path, cnn_model, mode='score', **kwargs):
    """
    Run gnina with a specific CNN model.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_path (str): Path to ligand file
        cnn_model (str): CNN model name
        mode (str): 'score' for scoring only, 'dock' for docking
        **kwargs: Additional parameters

    Returns:
        dict: Results including scores and timing
    """
    # Build command
    cmd = ['gnina', '-r', receptor_path, '-l', ligand_path, '--cnn', cnn_model]

    if mode == 'score':
        cmd.append('--score_only')
    else:
        # Add docking parameters
        if 'autobox_ligand' in kwargs:
            cmd.extend(['--autobox_ligand', kwargs['autobox_ligand']])
        elif 'center' in kwargs:
            cmd.extend(['--center_x', str(kwargs['center'][0])])
            cmd.extend(['--center_y', str(kwargs['center'][1])])
            cmd.extend(['--center_z', str(kwargs['center'][2])])
            cmd.extend(['--size_x', str(kwargs.get('size', [20, 20, 20])[0])])
            cmd.extend(['--size_y', str(kwargs.get('size', [20, 20, 20])[1])])
            cmd.extend(['--size_z', str(kwargs.get('size', [20, 20, 20])[2])])

        cmd.extend(['--num_modes', str(kwargs.get('num_modes', 3))])
        cmd.extend(['--exhaustiveness', str(kwargs.get('exhaustiveness', 8))])

    if kwargs.get('no_gpu', False):
        cmd.append('--no_gpu')

    try:
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        end_time = time.time()

        output = result.stdout + result.stderr
        parsed_results = parse_gnina_output_detailed(output)

        return {
            'model': cnn_model,
            'mode': mode,
            'success': True,
            'total_time': end_time - start_time,
            'ligand_name': Path(ligand_path).stem,
            **parsed_results
        }

    except subprocess.CalledProcessError as e:
        return {
            'model': cnn_model,
            'mode': mode,
            'success': False,
            'error': str(e),
            'ligand_name': Path(ligand_path).stem
        }
    except Exception as e:
        return {
            'model': cnn_model,
            'mode': mode,
            'success': False,
            'error': f"Unexpected error: {e}",
            'ligand_name': Path(ligand_path).stem
        }


def run_comprehensive_comparison(receptor_path, ligand_paths, cnn_models, mode='both', **kwargs):
    """
    Run comprehensive CNN model comparison.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_paths (list): List of ligand file paths
        cnn_models (list): List of CNN models to test
        mode (str): 'score', 'dock', or 'both'
        **kwargs: Additional parameters

    Returns:
        pandas.DataFrame: Comparison results
    """
    results = []

    modes_to_test = []
    if mode in ['score', 'both']:
        modes_to_test.append('score')
    if mode in ['dock', 'both']:
        modes_to_test.append('dock')

    total_runs = len(ligand_paths) * len(cnn_models) * len(modes_to_test)
    current_run = 0

    for ligand_path in ligand_paths:
        ligand_name = Path(ligand_path).stem
        print(f"\nProcessing ligand: {ligand_name}")

        for test_mode in modes_to_test:
            print(f"  Mode: {test_mode}")

            for cnn_model in cnn_models:
                current_run += 1
                print(f"    [{current_run}/{total_runs}] Testing {cnn_model}...", end=' ')

                result = run_cnn_comparison_single(
                    receptor_path, ligand_path, cnn_model, test_mode, **kwargs
                )

                if result['success']:
                    print(f"✓ Affinity: {result.get('affinity', 'N/A')}, "
                          f"CNN: {result.get('cnn_score', 'N/A')}, "
                          f"Time: {result['total_time']:.2f}s")
                else:
                    print(f"✗ Failed: {result.get('error', 'Unknown error')}")

                results.append(result)

    return pd.DataFrame(results)


def analyze_model_performance(results_df):
    """
    Analyze CNN model performance.

    Args:
        results_df (pandas.DataFrame): Results from model comparison

    Returns:
        dict: Performance analysis
    """
    successful = results_df[results_df['success'] == True].copy()

    if len(successful) == 0:
        return {'error': 'No successful results to analyze'}

    analysis = {}

    # Overall statistics
    analysis['total_runs'] = len(results_df)
    analysis['successful_runs'] = len(successful)
    analysis['success_rate'] = len(successful) / len(results_df) * 100

    # Performance by model
    model_stats = []
    for model in successful['model'].unique():
        model_data = successful[successful['model'] == model]

        stats = {
            'model': model,
            'runs': len(model_data),
            'success_rate': len(model_data) / len(results_df[results_df['model'] == model]) * 100,
            'mean_time': model_data['total_time'].mean(),
            'std_time': model_data['total_time'].std()
        }

        # Add score statistics if available
        if 'affinity' in model_data.columns and model_data['affinity'].notna().any():
            affinities = model_data['affinity'].dropna()
            stats.update({
                'mean_affinity': affinities.mean(),
                'std_affinity': affinities.std(),
                'best_affinity': affinities.min()
            })

        if 'cnn_score' in model_data.columns and model_data['cnn_score'].notna().any():
            cnn_scores = model_data['cnn_score'].dropna()
            stats.update({
                'mean_cnn_score': cnn_scores.mean(),
                'std_cnn_score': cnn_scores.std(),
                'best_cnn_score': cnn_scores.max()
            })

        model_stats.append(stats)

    analysis['model_performance'] = pd.DataFrame(model_stats)

    # Correlation analysis (if multiple ligands)
    if len(successful['ligand_name'].unique()) > 1:
        score_correlations = {}

        for model1 in successful['model'].unique():
            for model2 in successful['model'].unique():
                if model1 < model2:  # Avoid duplicates
                    data1 = successful[successful['model'] == model1]
                    data2 = successful[successful['model'] == model2]

                    # Merge on ligand_name and mode
                    merged = pd.merge(data1, data2, on=['ligand_name', 'mode'], suffixes=('_1', '_2'))

                    if len(merged) > 1 and 'affinity_1' in merged.columns:
                        corr = merged['affinity_1'].corr(merged['affinity_2'])
                        score_correlations[f"{model1}_vs_{model2}"] = corr

        analysis['score_correlations'] = score_correlations

    return analysis


def generate_comparison_plots(results_df, output_dir=None):
    """
    Generate comparison plots for CNN models.

    Args:
        results_df (pandas.DataFrame): Results dataframe
        output_dir (str): Directory to save plots (optional)

    Returns:
        list: List of generated plot files
    """
    successful = results_df[results_df['success'] == True].copy()

    if len(successful) == 0:
        print("No successful results to plot")
        return []

    plot_files = []

    # Set up plot style
    plt.style.use('default')
    sns.set_palette("husl")

    # 1. Runtime comparison
    if len(successful) > 0:
        plt.figure(figsize=(12, 6))

        # Box plot of runtimes by model
        plt.subplot(1, 2, 1)
        sns.boxplot(data=successful, x='model', y='total_time')
        plt.xticks(rotation=45, ha='right')
        plt.title('Runtime Comparison by CNN Model')
        plt.ylabel('Runtime (seconds)')

        # Bar plot of mean runtimes
        plt.subplot(1, 2, 2)
        runtime_means = successful.groupby('model')['total_time'].mean().sort_values()
        runtime_means.plot(kind='bar')
        plt.title('Mean Runtime by CNN Model')
        plt.ylabel('Mean Runtime (seconds)')
        plt.xticks(rotation=45, ha='right')

        plt.tight_layout()

        if output_dir:
            runtime_file = os.path.join(output_dir, 'cnn_runtime_comparison.png')
            plt.savefig(runtime_file, dpi=300, bbox_inches='tight')
            plot_files.append(runtime_file)
            print(f"Runtime plot saved: {runtime_file}")

        plt.show()

    # 2. Score comparison (if available)
    if 'affinity' in successful.columns and successful['affinity'].notna().any():
        plt.figure(figsize=(12, 6))

        # Affinity distribution by model
        plt.subplot(1, 2, 1)
        sns.violinplot(data=successful, x='model', y='affinity')
        plt.xticks(rotation=45, ha='right')
        plt.title('Affinity Distribution by CNN Model')
        plt.ylabel('Binding Affinity')

        # CNN Score vs Affinity correlation
        if 'cnn_score' in successful.columns and successful['cnn_score'].notna().any():
            plt.subplot(1, 2, 2)
            for model in successful['model'].unique():
                model_data = successful[successful['model'] == model]
                if len(model_data) > 1:
                    plt.scatter(model_data['affinity'], model_data['cnn_score'],
                              label=model, alpha=0.7)

            plt.xlabel('Binding Affinity')
            plt.ylabel('CNN Score')
            plt.title('CNN Score vs Affinity by Model')
            plt.legend()

        plt.tight_layout()

        if output_dir:
            score_file = os.path.join(output_dir, 'cnn_score_comparison.png')
            plt.savefig(score_file, dpi=300, bbox_inches='tight')
            plot_files.append(score_file)
            print(f"Score plot saved: {score_file}")

        plt.show()

    # 3. Performance summary heatmap
    if len(successful['model'].unique()) > 1:
        # Create performance matrix
        models = successful['model'].unique()
        metrics = []

        for model in models:
            model_data = successful[successful['model'] == model]
            metric_row = [
                model_data['total_time'].mean(),
                model_data['affinity'].mean() if 'affinity' in model_data.columns and model_data['affinity'].notna().any() else np.nan,
                model_data['cnn_score'].mean() if 'cnn_score' in model_data.columns and model_data['cnn_score'].notna().any() else np.nan
            ]
            metrics.append(metric_row)

        metrics_df = pd.DataFrame(metrics, index=models, columns=['Mean Runtime', 'Mean Affinity', 'Mean CNN Score'])

        # Normalize metrics for comparison
        metrics_normalized = metrics_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else x)

        plt.figure(figsize=(10, 8))
        sns.heatmap(metrics_normalized, annot=True, cmap='RdYlBu_r', center=0.5)
        plt.title('Normalized Performance Metrics by CNN Model')
        plt.tight_layout()

        if output_dir:
            heatmap_file = os.path.join(output_dir, 'cnn_performance_heatmap.png')
            plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plot_files.append(heatmap_file)
            print(f"Heatmap saved: {heatmap_file}")

        plt.show()

    return plot_files


def main():
    parser = argparse.ArgumentParser(
        description='Compare CNN models in gnina',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--receptor', '-r', required=True,
                       help='Receptor PDB file')
    parser.add_argument('--ligand', '-l',
                       help='Single ligand SDF/PDB file')
    parser.add_argument('--ligand_dir',
                       help='Directory containing multiple ligand files')
    parser.add_argument('--models', nargs='+',
                       default=['default', 'fast', 'dense'],
                       help='CNN models to compare')
    parser.add_argument('--mode', choices=['score', 'dock', 'both'],
                       default='score',
                       help='Comparison mode (default: score)')
    parser.add_argument('--autobox_ligand',
                       help='Reference ligand for autobox (required for docking mode)')
    parser.add_argument('--center', nargs=3, type=float,
                       help='Box center coordinates (x y z)')
    parser.add_argument('--size', nargs=3, type=float, default=[20, 20, 20],
                       help='Box size (default: 20 20 20)')
    parser.add_argument('--num_modes', type=int, default=3,
                       help='Number of docking modes (default: 3)')
    parser.add_argument('--output', '-o',
                       help='Output CSV file for results')
    parser.add_argument('--plot_dir',
                       help='Directory to save plots')
    parser.add_argument('--no_gpu', action='store_true',
                       help='Disable GPU acceleration')
    parser.add_argument('--list_models', action='store_true',
                       help='List available CNN models and exit')

    args = parser.parse_args()

    if args.list_models:
        models = get_available_cnn_models()
        print("Available CNN Models:")
        print("=" * 50)
        for model, description in models.items():
            print(f"{model:<25} : {description}")
        return

    # Validate inputs
    if not os.path.exists(args.receptor):
        print(f"Error: Receptor file {args.receptor} not found")
        sys.exit(1)

    # Collect ligand files
    ligand_paths = []
    if args.ligand:
        if not os.path.exists(args.ligand):
            print(f"Error: Ligand file {args.ligand} not found")
            sys.exit(1)
        ligand_paths = [args.ligand]
    elif args.ligand_dir:
        ligand_dir = Path(args.ligand_dir)
        ligand_paths = list(ligand_dir.glob('*.sdf')) + list(ligand_dir.glob('*.pdb'))
        if not ligand_paths:
            print(f"Error: No ligand files found in {args.ligand_dir}")
            sys.exit(1)
    else:
        print("Error: Must specify either --ligand or --ligand_dir")
        sys.exit(1)

    # Validate models
    available_models = get_available_cnn_models()
    invalid_models = [m for m in args.models if m not in available_models]
    if invalid_models:
        print(f"Error: Invalid models: {invalid_models}")
        print(f"Available models: {list(available_models.keys())}")
        sys.exit(1)

    # Setup parameters
    kwargs = {
        'no_gpu': args.no_gpu,
        'num_modes': args.num_modes
    }

    if args.mode in ['dock', 'both']:
        if args.autobox_ligand:
            if not os.path.exists(args.autobox_ligand):
                print(f"Error: Autobox ligand {args.autobox_ligand} not found")
                sys.exit(1)
            kwargs['autobox_ligand'] = args.autobox_ligand
        elif args.center:
            kwargs['center'] = args.center
            kwargs['size'] = args.size
        else:
            print("Error: For docking mode, must specify --autobox_ligand or --center")
            sys.exit(1)

    print(f"CNN Model Comparison")
    print(f"Receptor: {args.receptor}")
    print(f"Ligands: {len(ligand_paths)}")
    print(f"Models: {args.models}")
    print(f"Mode: {args.mode}")

    # Run comparison
    results_df = run_comprehensive_comparison(
        args.receptor,
        [str(p) for p in ligand_paths],
        args.models,
        args.mode,
        **kwargs
    )

    # Analyze results
    analysis = analyze_model_performance(results_df)

    if 'error' not in analysis:
        print("\n" + "="*60)
        print("CNN MODEL COMPARISON RESULTS")
        print("="*60)

        print(f"Total runs: {analysis['total_runs']}")
        print(f"Successful runs: {analysis['successful_runs']} ({analysis['success_rate']:.1f}%)")

        print("\nPERFORMANCE BY MODEL:")
        performance_df = analysis['model_performance']
        print(performance_df.to_string(index=False, float_format='%.3f'))

        if 'score_correlations' in analysis and analysis['score_correlations']:
            print("\nSCORE CORRELATIONS BETWEEN MODELS:")
            for pair, corr in analysis['score_correlations'].items():
                print(f"  {pair}: {corr:.3f}")

        # Save results
        if args.output:
            results_df.to_csv(args.output, index=False)
            print(f"\nDetailed results saved to: {args.output}")

        # Generate plots
        if args.plot_dir:
            os.makedirs(args.plot_dir, exist_ok=True)
            plot_files = generate_comparison_plots(results_df, args.plot_dir)
            print(f"\nPlots saved to: {args.plot_dir}")

    else:
        print(f"Analysis failed: {analysis['error']}")
        sys.exit(1)


if __name__ == '__main__':
    main()