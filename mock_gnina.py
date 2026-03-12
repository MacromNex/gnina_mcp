#!/usr/bin/env python3
"""
Mock gnina script for testing use cases when gnina binary is not available.
Simulates gnina output for different command line arguments.
"""
import sys
import random
import time

def main():
    args = sys.argv[1:]

    # Check for basic arguments
    if '--help' in args or '-h' in args:
        print("""gnina v1.3.2 (mock version for testing)
Usage: gnina [options]
  -r, --receptor arg       receptor (pdb or pdbqt)
  -l, --ligand arg         ligand (sdf, mol2, pdb, or pdbqt)
  --score_only             score only - do not dock
  --scoring arg            empirical scoring function
  --cnn arg                CNN model to use
  --autobox_ligand arg     automatically set box around this ligand
  --center_x arg           x coordinate of center
  --center_y arg           y coordinate of center
  --center_z arg           z coordinate of center
  --size_x arg             size in the x dimension (Angstroms)
  --size_y arg             size in the y dimension (Angstroms)
  --size_z arg             size in the z dimension (Angstroms)
  --out arg                output file path
  --num_modes arg          maximum number of modes to output (default 9)
  --exhaustiveness arg     exhaustiveness of the global search
  --flexres arg            flexible residues
  --flexdist arg           distance for flexible residues
  --cnn_scoring arg        CNN scoring mode (none/all/refinement)
""")
        return 0

    if '--version' in args:
        print("gnina v1.3.2 (mock version for testing)")
        return 0

    # Extract receptor and ligand from arguments
    receptor = None
    ligand = None
    score_only = False
    scoring_function = "default"
    cnn_model = "default"

    i = 0
    while i < len(args):
        if args[i] in ['-r', '--receptor'] and i + 1 < len(args):
            receptor = args[i + 1]
            i += 2
        elif args[i] in ['-l', '--ligand'] and i + 1 < len(args):
            ligand = args[i + 1]
            i += 2
        elif args[i] == '--score_only':
            score_only = True
            i += 1
        elif args[i] == '--scoring' and i + 1 < len(args):
            scoring_function = args[i + 1]
            i += 2
        elif args[i] == '--cnn' and i + 1 < len(args):
            cnn_model = args[i + 1]
            i += 2
        else:
            i += 1

    # Check for required arguments
    if not receptor:
        print("Error: receptor file required (-r/--receptor)", file=sys.stderr)
        return 1
    if not ligand:
        print("Error: ligand file required (-l/--ligand)", file=sys.stderr)
        return 1

    # Simulate processing time
    time.sleep(0.1)

    # Generate mock output
    if score_only:
        # Generate realistic mock scores
        affinity = round(random.uniform(-12.0, -6.0), 1)
        cnn_score = round(random.uniform(0.1, 0.95), 3)
        cnn_affinity = round(random.uniform(-10.0, -5.0), 1)

        print(f"""Receptor: {receptor}
Ligand: {ligand}
Scoring function: {scoring_function}
CNN model: {cnn_model}

Score only mode enabled.

Affinity: {affinity} (kcal/mol)
CNNscore: {cnn_score}
CNNaffinity: {cnn_affinity} (kcal/mol)

Pose 1:
   Affinity: {affinity}
   CNNscore: {cnn_score}
   CNNaffinity: {cnn_affinity}

Scoring completed successfully.""")
    else:
        # Mock docking output with multiple poses
        num_modes = 3  # Default number of poses

        # Try to extract num_modes from arguments
        for i in range(len(args)):
            if args[i] == '--num_modes' and i + 1 < len(args):
                try:
                    num_modes = int(args[i + 1])
                    num_modes = min(num_modes, 9)  # Cap at 9 for simplicity
                except ValueError:
                    pass
                break

        print(f"""Receptor: {receptor}
Ligand: {ligand}

Mode |   Affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------""")

        cnn_scores = []
        for i in range(num_modes):
            affinity = round(random.uniform(-11.0, -6.0), 1)
            rmsd_lb = round(random.uniform(0.0, 3.0), 3) if i > 0 else 0.0
            rmsd_ub = round(random.uniform(rmsd_lb, rmsd_lb + 2.0), 3) if i > 0 else 0.0
            cnn_score = round(random.uniform(0.3, 0.9), 3)
            cnn_affinity = round(random.uniform(-10.0, -5.0), 1)

            # Format compatible with the parser: mode affinity rmsd cnn_score
            print(f"   {i+1}       {affinity}      {rmsd_lb:.3f}     {cnn_score:.3f}")
            cnn_scores.append((cnn_score, cnn_affinity))

        print(f"""
CNN scoring:""")
        for i, (cnn_score, cnn_affinity) in enumerate(cnn_scores):
            print(f"Mode {i+1}: CNNscore = {cnn_score}, CNNaffinity = {cnn_affinity}")

        print(f"""
Docking completed successfully.""")

    return 0

if __name__ == "__main__":
    exit(main())