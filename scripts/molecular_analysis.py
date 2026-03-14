#!/usr/bin/env python3
"""
Script: molecular_analysis.py
Description: Molecular descriptor analysis and drug-likeness evaluation

Original Use Case: examples/use_case_6_python_api.py
Dependencies Removed: pygnina (optional, graceful fallback)

Usage:
    python scripts/molecular_analysis.py --ligand <ligand.sdf> --output <analysis.csv>

Example:
    python scripts/molecular_analysis.py --ligand examples/data/184l_lig.sdf --output results/analysis.csv
"""

import argparse
import sys
import os
import json
from pathlib import Path
from typing import Union, Optional, Dict, Any, List

import pandas as pd
import numpy as np

# Try to import RDKit for molecular manipulation
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Limited functionality.")

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
    "analysis": {
        "mode": "analyze",
        "include_descriptors": True,
        "drug_likeness": True
    },
    "descriptors": {
        "molecular_weight": True,
        "logp": True,
        "hbd_hba": True,
        "rotatable_bonds": True,
        "surface_area": True,
        "cyclic_peptide_specific": True
    },
    "drug_likeness": {
        "lipinski_rules": True,
        "veber_rules": True,
        "custom_thresholds": {
            "mw_max": 1000,
            "logp_max": 5,
            "hbd_max": 10,
            "hba_max": 10
        }
    },
    "output": {
        "format": "csv",
        "include_report": True,
        "detailed_analysis": True
    }
}

def calculate_molecular_descriptors(mol: 'Chem.Mol') -> Dict[str, Any]:
    """Calculate molecular descriptors using RDKit."""
    if not RDKIT_AVAILABLE:
        return {"error": "RDKit not available"}

    descriptors = {}

    try:
        # Basic descriptors
        descriptors['molecular_weight'] = Descriptors.MolWt(mol)
        descriptors['logp'] = Descriptors.MolLogP(mol)
        descriptors['hbd'] = Descriptors.NumHDonors(mol)
        descriptors['hba'] = Descriptors.NumHAcceptors(mol)
        descriptors['rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
        descriptors['rings'] = rdMolDescriptors.CalcNumRings(mol)
        descriptors['aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol)
        descriptors['heavy_atoms'] = mol.GetNumHeavyAtoms()

        # Surface area (approximation)
        try:
            descriptors['tpsa'] = Descriptors.TPSA(mol)
        except:
            descriptors['tpsa'] = None

        # Complexity measures
        descriptors['complexity'] = Descriptors.BertzCT(mol)

        # Cyclic peptide specific
        ring_info = mol.GetRingInfo()
        descriptors['is_cyclic'] = ring_info.NumRings() > 0
        descriptors['largest_ring'] = max([len(ring) for ring in ring_info.AtomRings()]) if ring_info.AtomRings() else 0

    except Exception as e:
        descriptors['error'] = str(e)

    return descriptors

def evaluate_drug_likeness(descriptors: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
    """Evaluate drug-likeness using various rules."""
    evaluation = {}

    # Lipinski Rule of Five
    if config['drug_likeness']['lipinski_rules']:
        lipinski_violations = 0
        lipinski_details = {}

        if descriptors.get('molecular_weight', 0) > 500:
            lipinski_violations += 1
            lipinski_details['mw_violation'] = True
        if descriptors.get('logp', 0) > 5:
            lipinski_violations += 1
            lipinski_details['logp_violation'] = True
        if descriptors.get('hbd', 0) > 5:
            lipinski_violations += 1
            lipinski_details['hbd_violation'] = True
        if descriptors.get('hba', 0) > 10:
            lipinski_violations += 1
            lipinski_details['hba_violation'] = True

        evaluation['lipinski_violations'] = lipinski_violations
        evaluation['lipinski_compliant'] = lipinski_violations <= 1  # Allow 1 violation
        evaluation['lipinski_details'] = lipinski_details

    # Veber Rules
    if config['drug_likeness']['veber_rules']:
        veber_violations = 0
        veber_details = {}

        if descriptors.get('rotatable_bonds', 0) > 10:
            veber_violations += 1
            veber_details['rotbond_violation'] = True
        if descriptors.get('tpsa', 0) > 140:
            veber_violations += 1
            veber_details['tpsa_violation'] = True

        evaluation['veber_violations'] = veber_violations
        evaluation['veber_compliant'] = veber_violations == 0
        evaluation['veber_details'] = veber_details

    # Custom thresholds for cyclic peptides
    custom = config['drug_likeness']['custom_thresholds']
    custom_violations = 0

    if descriptors.get('molecular_weight', 0) > custom['mw_max']:
        custom_violations += 1
    if descriptors.get('logp', 0) > custom['logp_max']:
        custom_violations += 1
    if descriptors.get('hbd', 0) > custom['hbd_max']:
        custom_violations += 1
    if descriptors.get('hba', 0) > custom['hba_max']:
        custom_violations += 1

    evaluation['custom_violations'] = custom_violations
    evaluation['custom_compliant'] = custom_violations <= 2  # More lenient for cyclic peptides

    return evaluation

def run_molecular_analysis(
    ligand_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform molecular descriptor analysis and drug-likeness evaluation.

    Args:
        ligand_file: Path to ligand file (SDF, MOL, etc.)
        output_file: Path to save analysis results
        config: Configuration dict
        **kwargs: Override specific config parameters

    Returns:
        Dict containing analysis results
    """
    # Setup
    ligand_file = Path(ligand_file)
    config = _deep_merge(DEFAULT_CONFIG, config or {})
    if kwargs:
        config = _deep_merge(config, kwargs)

    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for molecular analysis")

    if not ligand_file.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")

    # Load molecules
    if ligand_file.suffix.lower() == '.sdf':
        suppl = Chem.SDMolSupplier(str(ligand_file))
        molecules = [mol for mol in suppl if mol is not None]
    elif ligand_file.suffix.lower() in ['.mol', '.mol2']:
        mol = Chem.MolFromMolFile(str(ligand_file))
        molecules = [mol] if mol is not None else []
    else:
        raise ValueError(f"Unsupported file format: {ligand_file.suffix}")

    if not molecules:
        raise ValueError("No valid molecules found in file")

    # Analyze each molecule
    results = []
    for i, mol in enumerate(molecules):
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i+1}"

        # Calculate descriptors
        descriptors = calculate_molecular_descriptors(mol)

        # Evaluate drug-likeness
        drug_likeness = evaluate_drug_likeness(descriptors, config)

        # Combine results
        result = {
            "molecule": mol_name,
            "file": ligand_file.name,
            **descriptors,
            **drug_likeness
        }

        results.append(result)

    # Generate summary statistics
    if results:
        summary = {
            "total_molecules": len(results),
            "lipinski_compliant": sum(1 for r in results if r.get('lipinski_compliant', False)),
            "veber_compliant": sum(1 for r in results if r.get('veber_compliant', False)),
            "custom_compliant": sum(1 for r in results if r.get('custom_compliant', False)),
            "cyclic_molecules": sum(1 for r in results if r.get('is_cyclic', False))
        }

        # Calculate mean descriptors
        numeric_descriptors = ['molecular_weight', 'logp', 'hbd', 'hba', 'rotatable_bonds', 'tpsa']
        for desc in numeric_descriptors:
            values = [r[desc] for r in results if r.get(desc) is not None]
            if values:
                summary[f'mean_{desc}'] = np.mean(values)
                summary[f'std_{desc}'] = np.std(values)
    else:
        summary = {"error": "No results to analyze"}

    # Save results
    output_path = None
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save detailed results
        df = pd.DataFrame(results)
        df.to_csv(output_path, index=False, float_format='%.3f')

        # Save report
        if config['output']['include_report']:
            report_path = output_path.parent / f"{output_path.stem}_report.txt"
            with open(report_path, 'w') as f:
                f.write("="*60 + "\n")
                f.write("MOLECULAR ANALYSIS REPORT\n")
                f.write("="*60 + "\n")
                for key, value in summary.items():
                    f.write(f"{key}: {value}\n")

    return {
        "results": results,
        "summary": summary,
        "output_file": str(output_path) if output_path else None,
        "metadata": {
            "ligand_file": str(ligand_file),
            "config": config
        }
    }

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ligand', '-l', required=True, help='Ligand file (SDF, MOL)')
    parser.add_argument('--output', '-o', help='Output CSV file')
    parser.add_argument('--config', '-c', help='Config file (JSON)')

    args = parser.parse_args()

    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    try:
        result = run_molecular_analysis(
            ligand_file=args.ligand,
            output_file=args.output,
            config=config
        )

        print("="*60)
        print("MOLECULAR ANALYSIS RESULTS")
        print("="*60)

        # Print summary
        summary = result['summary']
        for key, value in summary.items():
            print(f"{key}: {value}")

        # Print detailed results
        if result['results']:
            print("\nDETAILED RESULTS:")
            df = pd.DataFrame(result['results'])
            key_columns = ['molecule', 'molecular_weight', 'logp', 'hbd', 'hba', 'lipinski_compliant']
            display_columns = [col for col in key_columns if col in df.columns]
            print(df[display_columns].to_string(index=False, float_format='%.2f'))

        if result['output_file']:
            print(f"\nDetailed results saved to: {result['output_file']}")

        return result

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()