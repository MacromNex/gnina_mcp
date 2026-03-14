# Example usages:
System 1: 3rod (receptor + ligand + reference)

  # Score existing pose
  gnina -r examples/data/3rod_rec.pdb -l examples/data/3rod_lig.pdb --score_only

  # Dock using reference ligand to define binding site
  gnina -r examples/data/3rod_rec.pdb -l examples/data/3rod_lig.pdb --autobox_ligand examples/data/3rod_rec_ref.pdb -o results/3rod_docked.sdf.gz

  # Minimize existing pose
  gnina -r examples/data/3rod_rec.pdb -l examples/data/3rod_lig.pdb --minimize -o results/3rod_minimized.sdf.gz

  System 2: 10gs (receptor + ligand + reference)

  # Score existing pose
  gnina -r examples/data/10gs_rec.pdb -l examples/data/10gs_lig.sdf --score_only

  # Dock using reference ligand to define binding site
  gnina -r examples/data/10gs_rec.pdb -l examples/data/10gs_lig.sdf --autobox_ligand examples/data/10gs_rec_ref.pdb -o results/10gs_docked.sdf.gz

  # Dock with Vinardo scoring (no CNN)
  gnina -r examples/data/10gs_rec.pdb -l examples/data/10gs_lig.sdf --autobox_ligand examples/data/10gs_rec_ref.pdb --scoring vinardo --cnn_scoring none -o
  results/10gs_vinardo.sdf.gz

  System 3: 184l (receptor + ligand + 3 chain references)

  # Score existing pose
  gnina -r examples/data/184l_rec.pdb -l examples/data/184l_lig.sdf --score_only

  # Dock using reference ligand to define binding site
  gnina -r examples/data/184l_rec.pdb -l examples/data/184l_lig.sdf --autobox_ligand examples/data/184l_rec_ref.pdb -o results/184l_docked.sdf.gz

  # Dock with CNN refinement (slower but more accurate)
  gnina -r examples/data/184l_rec.pdb -l examples/data/184l_lig.sdf --autobox_ligand examples/data/184l_rec_ref.pdb --cnn_scoring refinement -o
  results/184l_cnn_refined.sdf.gz

  # Flexible docking using reference to define flexible residues
  gnina -r examples/data/184l_rec.pdb -l examples/data/184l_lig.sdf --autobox_ligand examples/data/184l_rec_ref.pdb --flexdist_ligand examples/data/184l_rec_ref.pdb
  --flexdist 3.5 -o results/184l_flex_docked.sdf.gz


