# Changelog

## 1.0.0 - 2026-11-24

- Replaced the nested `data/data/` layout. Results are read from `data/docking_results.csv`, and the analysis and plotting scripts resolve paths from the repository root.
- Added a scripted AutoDock Vina 1.2.7 workflow that downloads PDB 6LU7 and the five PubChem structures by accession, prepares the receptor and ligands, and docks them in one fixed box (30 Å, seed 42, exhaustiveness 16).
- Replaced the previous score table with the scores from that run. Kaempferol −7.745, berberine −7.679, quercetin −7.518, curcumin −6.910, andrographolide −6.883 kcal/mol.
- Redocked the co-crystallised inhibitor (chain C of 6LU7, covalent C–S bond omitted). Symmetry-corrected heavy-atom RMSD in the crystal frame is 7.266 Å. The crystal coordinates score −3.470 kcal/mol; the top search pose scores −8.520 kcal/mol.
- Documented that these Vina scores are not measured binding affinities.
- Added PDB 6W63 / X77 as the non-covalent redocking control (symmetry-corrected heavy-atom RMSD 0.915 Å, which meets 2 Å). Kept the 6LU7 / N3 result (7.266 Å) as a covalent-complex negative control.
- Docked the same five phytochemicals into 6W63 and wrote a second score table. The rank order is not the same as on 6LU7.
- Added a pinned conda environment, pytest coverage of preparation and parsing, and a lightweight GitHub Actions workflow that does not launch Vina.
