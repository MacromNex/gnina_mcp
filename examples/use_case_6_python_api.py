#!/usr/bin/env python3
"""
Use Case 6: Python API Usage with PyGnina
==========================================

This script demonstrates how to use the pygnina Python interface for integrating
gnina functionality directly into Python workflows. This is useful for building
automated pipelines, custom analysis tools, and integration with other Python
molecular modeling libraries.

Use case: Programmatic access to gnina for custom workflows, automated screening
pipelines, and integration with machine learning or optimization frameworks.

Requirements:
- pygnina Python package (built with gnina)
- RDKit for molecular manipulation
- Receptor PDB file(s)
- Ligand SDF file(s)

Example Usage:
    python use_case_6_python_api.py --receptor data/184l_rec.pdb --ligand data/184l_lig.sdf
    python use_case_6_python_api.py --receptor data/3rod_rec.pdb --ligands data/3rod_lig.pdb data/10gs_lig.sdf --mode dock
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
import time
from typing import List, Dict, Any, Tuple, Optional

# Check for pygnina availability
try:
    import pygnina
    PYGNINA_AVAILABLE = True
    print("pygnina module loaded successfully")
except ImportError:
    PYGNINA_AVAILABLE = False
    print("Warning: pygnina module not available. This may be expected if gnina is not built with Python support.")

# Try to import RDKit for molecular manipulation
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Some functionality will be limited.")


class GninaInterface:
    """
    Python interface wrapper for gnina functionality.
    """

    def __init__(self):
        """Initialize the gnina interface."""
        self.gnina_instance = None
        self.receptor_loaded = False
        self.receptor_path = None

        if PYGNINA_AVAILABLE:
            try:
                self.gnina_instance = pygnina.GNINA()
                print("PyGnina instance created successfully")
            except Exception as e:
                print(f"Error creating pygnina instance: {e}")
                self.gnina_instance = None

    def load_receptor(self, receptor_path: str, format_type: str = 'pdb') -> bool:
        """
        Load a receptor structure.

        Args:
            receptor_path (str): Path to receptor file
            format_type (str): File format ('pdb', 'mol2', etc.)

        Returns:
            bool: Success status
        """
        if not self.gnina_instance:
            print("Error: PyGnina not available")
            return False

        if not os.path.exists(receptor_path):
            print(f"Error: Receptor file not found: {receptor_path}")
            return False

        try:
            # Read receptor file content
            with open(receptor_path, 'r') as f:
                receptor_content = f.read()

            # Set receptor in pygnina
            self.gnina_instance.set_receptor(receptor_content, format_type)
            self.receptor_loaded = True
            self.receptor_path = receptor_path

            print(f"Receptor loaded: {Path(receptor_path).name}")
            return True

        except Exception as e:
            print(f"Error loading receptor: {e}")
            return False

    def score_ligand(self, ligand_path: str, format_type: str = 'sdf') -> Optional[Dict[str, float]]:
        """
        Score a single ligand against the loaded receptor.

        Args:
            ligand_path (str): Path to ligand file
            format_type (str): File format ('sdf', 'pdb', 'mol2')

        Returns:
            dict: Scoring results or None if failed
        """
        if not self.receptor_loaded:
            print("Error: No receptor loaded")
            return None

        if not os.path.exists(ligand_path):
            print(f"Error: Ligand file not found: {ligand_path}")
            return None

        try:
            # Read ligand file content
            with open(ligand_path, 'r') as f:
                ligand_content = f.read()

            # Score ligand
            start_time = time.time()
            result = self.gnina_instance.score(ligand_content, format_type)
            end_time = time.time()

            # Extract scores
            scores = {
                'ligand_name': Path(ligand_path).stem,
                'energy': result.energy(),
                'cnn_score': result.cnnscore(),
                'cnn_affinity': result.cnnaffinity(),
                'scoring_time': end_time - start_time
            }

            return scores

        except Exception as e:
            print(f"Error scoring ligand {Path(ligand_path).name}: {e}")
            return None

    def score_multiple_ligands(self, ligand_paths: List[str]) -> pd.DataFrame:
        """
        Score multiple ligands.

        Args:
            ligand_paths (list): List of ligand file paths

        Returns:
            pandas.DataFrame: Scoring results
        """
        if not self.receptor_loaded:
            print("Error: No receptor loaded")
            return pd.DataFrame()

        results = []
        total = len(ligand_paths)

        print(f"Scoring {total} ligands against {Path(self.receptor_path).name}")

        for i, ligand_path in enumerate(ligand_paths, 1):
            ligand_name = Path(ligand_path).name

            print(f"[{i}/{total}] Scoring {ligand_name}...", end=' ')

            scores = self.score_ligand(ligand_path)

            if scores:
                results.append(scores)
                print(f"✓ Energy: {scores['energy']:.3f}, CNN: {scores['cnn_score']:.3f}")
            else:
                print("✗ Failed")

        return pd.DataFrame(results)


class MolecularAnalyzer:
    """
    Analyzer for molecular properties and descriptors.
    """

    @staticmethod
    def calculate_descriptors(ligand_path: str) -> Optional[Dict[str, float]]:
        """
        Calculate molecular descriptors for a ligand.

        Args:
            ligand_path (str): Path to ligand SDF file

        Returns:
            dict: Molecular descriptors or None if failed
        """
        if not RDKIT_AVAILABLE:
            return None

        try:
            # Read molecule
            mol = Chem.SDMolSupplier(ligand_path)[0]

            if mol is None:
                return None

            # Calculate descriptors
            descriptors = {
                'molecular_weight': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
                'tpsa': Descriptors.TPSA(mol),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'formal_charge': Chem.rdmolops.GetFormalCharge(mol)
            }

            # Calculate cyclic peptide specific descriptors
            if mol.GetNumHeavyAtoms() > 20:  # Likely peptide
                # Count amino acid residues (simplified)
                pattern = Chem.MolFromSmarts('[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]')
                if pattern:
                    aa_matches = mol.GetSubstructMatches(pattern)
                    descriptors['peptide_residues'] = len(aa_matches)

                # Check for cyclicity
                ring_info = mol.GetRingInfo()
                descriptors['ring_count'] = ring_info.NumRings()
                descriptors['largest_ring'] = max([len(ring) for ring in ring_info.AtomRings()]) if ring_info.AtomRings() else 0

            return descriptors

        except Exception as e:
            print(f"Error calculating descriptors for {Path(ligand_path).name}: {e}")
            return None

    @staticmethod
    def analyze_drug_likeness(descriptors: Dict[str, float]) -> Dict[str, Any]:
        """
        Analyze drug-likeness using Lipinski's Rule of Five and other criteria.

        Args:
            descriptors (dict): Molecular descriptors

        Returns:
            dict: Drug-likeness analysis
        """
        analysis = {}

        # Lipinski's Rule of Five
        lipinski_violations = 0
        lipinski_criteria = {}

        # Molecular weight <= 500 Da
        mw_ok = descriptors.get('molecular_weight', 0) <= 500
        lipinski_criteria['molecular_weight'] = mw_ok
        if not mw_ok:
            lipinski_violations += 1

        # LogP <= 5
        logp_ok = descriptors.get('logp', 0) <= 5
        lipinski_criteria['logp'] = logp_ok
        if not logp_ok:
            lipinski_violations += 1

        # HBD <= 5
        hbd_ok = descriptors.get('hbd', 0) <= 5
        lipinski_criteria['hbd'] = hbd_ok
        if not hbd_ok:
            lipinski_violations += 1

        # HBA <= 10
        hba_ok = descriptors.get('hba', 0) <= 10
        lipinski_criteria['hba'] = hba_ok
        if not hba_ok:
            lipinski_violations += 1

        analysis['lipinski'] = {
            'criteria': lipinski_criteria,
            'violations': lipinski_violations,
            'drug_like': lipinski_violations <= 1
        }

        # Veber criteria (oral bioavailability)
        veber_ok = (descriptors.get('tpsa', 0) <= 140 and
                   descriptors.get('rotatable_bonds', 0) <= 10)

        analysis['veber'] = {
            'oral_bioavailable': veber_ok,
            'tpsa_ok': descriptors.get('tpsa', 0) <= 140,
            'rotbonds_ok': descriptors.get('rotatable_bonds', 0) <= 10
        }

        # Cyclic peptide specific analysis
        if descriptors.get('peptide_residues', 0) > 0:
            analysis['peptide'] = {
                'is_peptide': True,
                'residue_count': descriptors.get('peptide_residues', 0),
                'is_cyclic': descriptors.get('largest_ring', 0) > 6,
                'largest_ring': descriptors.get('largest_ring', 0)
            }

            # Cyclic peptide drug-likeness (Beyond Rule of Five)
            analysis['peptide']['bro5_compliant'] = (
                descriptors.get('molecular_weight', 0) <= 2000 and
                descriptors.get('hbd', 0) <= 12 and
                descriptors.get('hba', 0) <= 18
            )

        return analysis


def run_integrated_analysis(receptor_path: str, ligand_paths: List[str]) -> pd.DataFrame:
    """
    Run integrated analysis combining gnina scoring and molecular descriptors.

    Args:
        receptor_path (str): Path to receptor PDB file
        ligand_paths (list): List of ligand file paths

    Returns:
        pandas.DataFrame: Combined analysis results
    """
    print("Running integrated molecular analysis...")

    # Initialize gnina interface
    gnina = GninaInterface()

    if not gnina.load_receptor(receptor_path):
        return pd.DataFrame()

    # Score ligands
    scoring_results = gnina.score_multiple_ligands(ligand_paths)

    if scoring_results.empty:
        print("No scoring results obtained")
        return pd.DataFrame()

    # Calculate molecular descriptors
    descriptor_results = []
    analyzer = MolecularAnalyzer()

    print("\nCalculating molecular descriptors...")

    for ligand_path in ligand_paths:
        ligand_name = Path(ligand_path).stem

        descriptors = analyzer.calculate_descriptors(ligand_path)

        if descriptors:
            descriptors['ligand_name'] = ligand_name

            # Analyze drug-likeness
            drug_analysis = analyzer.analyze_drug_likeness(descriptors)

            # Flatten drug-likeness results
            descriptors['lipinski_violations'] = drug_analysis['lipinski']['violations']
            descriptors['drug_like'] = drug_analysis['lipinski']['drug_like']
            descriptors['oral_bioavailable'] = drug_analysis['veber']['oral_bioavailable']

            if 'peptide' in drug_analysis:
                descriptors['is_peptide'] = drug_analysis['peptide']['is_peptide']
                descriptors['is_cyclic'] = drug_analysis['peptide']['is_cyclic']
                descriptors['bro5_compliant'] = drug_analysis['peptide']['bro5_compliant']

            descriptor_results.append(descriptors)

    descriptor_df = pd.DataFrame(descriptor_results)

    # Merge scoring and descriptor results
    if not descriptor_df.empty:
        combined_results = pd.merge(scoring_results, descriptor_df, on='ligand_name', how='outer')
    else:
        combined_results = scoring_results

    return combined_results


def generate_analysis_report(results_df: pd.DataFrame, output_file: Optional[str] = None):
    """
    Generate comprehensive analysis report.

    Args:
        results_df (pandas.DataFrame): Analysis results
        output_file (str): Output file path (optional)
    """
    if results_df.empty:
        print("No results to analyze")
        return

    report_lines = []
    report_lines.append("INTEGRATED MOLECULAR ANALYSIS REPORT")
    report_lines.append("=" * 60)

    # Summary statistics
    report_lines.append(f"\nANALYZED COMPOUNDS: {len(results_df)}")

    # Scoring statistics
    if 'energy' in results_df.columns:
        energies = results_df['energy'].dropna()
        if not energies.empty:
            report_lines.append(f"\nSCORING STATISTICS:")
            report_lines.append(f"  Best binding energy: {energies.min():.3f}")
            report_lines.append(f"  Mean binding energy: {energies.mean():.3f}")
            report_lines.append(f"  Energy range: {energies.max() - energies.min():.3f}")

    if 'cnn_score' in results_df.columns:
        cnn_scores = results_df['cnn_score'].dropna()
        if not cnn_scores.empty:
            report_lines.append(f"  Best CNN score: {cnn_scores.max():.3f}")
            report_lines.append(f"  Mean CNN score: {cnn_scores.mean():.3f}")

    # Drug-likeness statistics
    if 'drug_like' in results_df.columns:
        drug_like_count = results_df['drug_like'].sum()
        report_lines.append(f"\nDRUG-LIKENESS ANALYSIS:")
        report_lines.append(f"  Lipinski compliant: {drug_like_count}/{len(results_df)} ({drug_like_count/len(results_df)*100:.1f}%)")

        if 'oral_bioavailable' in results_df.columns:
            oral_count = results_df['oral_bioavailable'].sum()
            report_lines.append(f"  Oral bioavailable: {oral_count}/{len(results_df)} ({oral_count/len(results_df)*100:.1f}%)")

    # Peptide analysis
    if 'is_peptide' in results_df.columns:
        peptides = results_df[results_df['is_peptide'] == True]
        if not peptides.empty:
            report_lines.append(f"\nPEPTIDE ANALYSIS:")
            report_lines.append(f"  Peptides identified: {len(peptides)}")

            if 'is_cyclic' in peptides.columns:
                cyclic_count = peptides['is_cyclic'].sum()
                report_lines.append(f"  Cyclic peptides: {cyclic_count}")

            if 'bro5_compliant' in peptides.columns:
                bro5_count = peptides['bro5_compliant'].sum()
                report_lines.append(f"  Beyond RO5 compliant: {bro5_count}")

    # Top compounds
    if 'energy' in results_df.columns:
        top_compounds = results_df.nsmallest(5, 'energy')[['ligand_name', 'energy', 'cnn_score']].round(3)
        report_lines.append(f"\nTOP 5 COMPOUNDS BY BINDING ENERGY:")
        report_lines.append(top_compounds.to_string(index=False))

    report_text = "\n".join(report_lines)
    print(report_text)

    # Save to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
        print(f"\nReport saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Python API usage example for gnina integration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--receptor', '-r', required=True,
                       help='Receptor PDB file')
    parser.add_argument('--ligand', '-l',
                       help='Single ligand file')
    parser.add_argument('--ligands', nargs='+',
                       help='Multiple ligand files')
    parser.add_argument('--ligand_dir',
                       help='Directory containing ligand files')
    parser.add_argument('--mode', choices=['score', 'analyze'],
                       default='analyze',
                       help='Analysis mode (default: analyze)')
    parser.add_argument('--output', '-o',
                       help='Output CSV file for results')
    parser.add_argument('--report',
                       help='Output file for analysis report')
    parser.add_argument('--check_dependencies', action='store_true',
                       help='Check available dependencies and exit')

    args = parser.parse_args()

    if args.check_dependencies:
        print("Checking dependencies:")
        print(f"  PyGnina available: {PYGNINA_AVAILABLE}")
        print(f"  RDKit available: {RDKIT_AVAILABLE}")

        if not PYGNINA_AVAILABLE:
            print("\nTo enable PyGnina:")
            print("  1. Build gnina with Python support")
            print("  2. Install the pygnina package")

        if not RDKIT_AVAILABLE:
            print("\nTo enable RDKit:")
            print("  conda install -c conda-forge rdkit")

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
    elif args.ligands:
        for ligand in args.ligands:
            if not os.path.exists(ligand):
                print(f"Error: Ligand file {ligand} not found")
                sys.exit(1)
        ligand_paths = args.ligands
    elif args.ligand_dir:
        ligand_dir = Path(args.ligand_dir)
        ligand_paths = list(ligand_dir.glob('*.sdf')) + list(ligand_dir.glob('*.pdb'))
        ligand_paths = [str(p) for p in ligand_paths]
    else:
        print("Error: Must specify ligands using --ligand, --ligands, or --ligand_dir")
        sys.exit(1)

    if not ligand_paths:
        print("Error: No ligand files found")
        sys.exit(1)

    print(f"Python API Analysis Mode: {args.mode}")
    print(f"Receptor: {args.receptor}")
    print(f"Ligands: {len(ligand_paths)} files")

    if not PYGNINA_AVAILABLE:
        print("\nWarning: PyGnina not available. Falling back to descriptor analysis only.")

        if args.mode == 'score':
            print("Error: Scoring mode requires PyGnina")
            sys.exit(1)

        # Run descriptor-only analysis
        analyzer = MolecularAnalyzer()
        results = []

        for ligand_path in ligand_paths:
            descriptors = analyzer.calculate_descriptors(ligand_path)
            if descriptors:
                drug_analysis = analyzer.analyze_drug_likeness(descriptors)
                descriptors.update({
                    'lipinski_violations': drug_analysis['lipinski']['violations'],
                    'drug_like': drug_analysis['lipinski']['drug_like'],
                    'oral_bioavailable': drug_analysis['veber']['oral_bioavailable']
                })
                results.append(descriptors)

        results_df = pd.DataFrame(results)

    else:
        # Run full analysis
        if args.mode == 'score':
            # Score-only mode
            gnina = GninaInterface()
            if gnina.load_receptor(args.receptor):
                results_df = gnina.score_multiple_ligands(ligand_paths)
            else:
                print("Failed to load receptor")
                sys.exit(1)

        else:
            # Full analysis mode
            results_df = run_integrated_analysis(args.receptor, ligand_paths)

    # Display and save results
    if not results_df.empty:
        print(f"\nAnalysis completed for {len(results_df)} compounds")

        # Save detailed results
        if args.output:
            results_df.to_csv(args.output, index=False)
            print(f"Detailed results saved to: {args.output}")

        # Generate report
        generate_analysis_report(results_df, args.report)

    else:
        print("No results obtained")
        sys.exit(1)


if __name__ == '__main__':
    main()