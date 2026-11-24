# Protocol

Virtual screening of five plant-derived compounds against SARS-CoV-2 main protease, with a redocking control. Numerical results are the run archived in this repository (2026-11-24). They are not copied from an earlier table.

## Structures

1. Download PDB 6LU7 from `https://files.rcsb.org/download/6LU7.pdb`. The file used for the published run is `structures/raw/6LU7.pdb` (SHA-256 in `results/run_metadata.json`). Resolution 2.16 Å, one chain of the protease (A, residues 1–306) and the inhibitor as chain C.
2. Receptor: ATOM records of chain A only, altloc blank or A, no hydrogens from the deposit. `mk_prepare_receptor.py` (Meeko 0.8.0) writes `structures/receptor/6LU7_receptor.pdbqt` from its residue templates.
3. Inhibitor: heavy atoms of chain C. Intra-residue bonds of 02J, PJE and 010 come from `data/ccd/*.cif`. ALA, VAL and LEU use standard amino-acid bonds. Peptide bonds join consecutive residues. LINK records that stay inside chain C are added, including the ester from PJE C to 010 O. The LINK from Cys145 SG to PJE C20 is not added. Hydrogens are added with RDKit and relaxed with UFF while heavy atoms are held at the crystal coordinates.
4. Phytochemicals: PubChem 3D SDF for the CID in `data/ligands.csv` (`record_type=3d`). If a 3D record were missing, the script would embed a 2D record with ETKDGv3 and seed 42. All five published ligands used the 3D record. Meeko writes the PDBQT and a `REMARK INDEX MAP` so pose coordinates can be put back onto the RDKit atom order. Nonpolar hydrogens are united atoms in the PDBQT; the RMSD uses heavy atoms only.

## Docking

AutoDock Vina 1.2.7, scoring function `vina`, rigid receptor.

The config `docking/vina_config.txt` is written by the script, not by hand:

- centre = centroid of the 49 crystal heavy atoms, rounded to 3 decimals: −10.729, 12.418, 68.816
- cubic edge = smallest even integer that is ≥ 22 Å and ≥ 2 × (maximum absolute deviation from the centroid + 4 Å) = 30 Å
- exhaustiveness 16, nine modes, energy range 3, mode separation 1 Å, seed 42, 4 threads

The crystal N3 PDBQT is scored with `--score_only`, then docked. Each phytochemical is docked with the same config. The best pose is model 1. The score in `data/docking_results.csv` is the `REMARK VINA RESULT` value from that model, checked against the log.

## RMSD

`src/mpro_docking/rmsd.py` computes the heavy-atom RMSD in the shared frame. It does not rotate or translate one pose onto the other. Symmetry-corrected RMSD is the minimum over graph automorphisms; direct RMSD keeps the atom order. The pose-reproduction cutoff used in this repository is 2 Å.

**PDB 6W63 / X77 (the control).** Non-covalent: no LINK record, closest crystal contact 2.796 Å. Box from the same rule: 22 Å cube centred at (−20.096, 18.844, −27.394). Symmetry-corrected RMSD 0.915 Å, direct RMSD 1.168 Å. Both meet 2 Å. Crystal pose score −8.076 kcal/mol; top search pose −8.507 kcal/mol.

**PDB 6LU7 / N3 (negative control).** Covalent complex. Symmetry-corrected RMSD 7.266 Å, direct RMSD 7.295 Å. The crystal pose scores −3.470 kcal/mol and the top search pose scores −8.520 kcal/mol. PJE C20 moves from 1.793 Å to 5.121 Å from Cys145 SG. Standard Vina has no term for that bond, so this failure is expected and is not used as evidence about the phytochemical poses.

## Two screening tables

The five PubChem ligands are docked into both receptors with seed 42 and exhaustiveness 16. `data/docking_results_6w63.csv` is the screen on the receptor whose co-crystal pose was recovered. `data/docking_results.csv` is the screen on 6LU7. The rank order is not the same. Neither table is a binding affinity.

## Interactions

From the docked PDBQT coordinates and the receptor PDBQT:

- hydrogen bond: donor and acceptor are N or O, H···A ≤ 2.5 Å, D–H···A ≥ 120°
- contact residue: any protein heavy atom within 4.0 Å of a ligand heavy atom

## What is not claimed

Vina scores are not binding affinities. The log line "Estimated Free Energy of Binding" is the program's name for its score. Ranks are only for this receptor, this box, this seed and these input conformers.
