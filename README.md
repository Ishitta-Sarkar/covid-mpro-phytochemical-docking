# COVID-19 Mpro Phytochemical Docking

A computational molecular docking workflow developed to explore the virtual screening of selected plant-derived compounds against the SARS-CoV-2 main protease.

## Project Overview

SARS-CoV-2 main protease is an important viral enzyme involved in polyprotein processing and viral replication. Because of its functional importance, it has been widely studied as a potential target in computational drug discovery and virtual screening studies.

This project presents a structured molecular docking workflow for selected phytochemicals against SARS-CoV-2 main protease using publicly available structural biology resources and docking-based analysis.

## Aim

To evaluate selected plant-derived compounds against SARS-CoV-2 main protease using a computational molecular docking workflow.

## Objectives

- Identify SARS-CoV-2 main protease as a target protein.
- Select plant-derived compounds reported in biomedical literature.
- Organize ligand and docking result datasets.
- Rank ligands based on predicted binding affinity.
- Summarize the virtual screening results in a reproducible format.

## Repository Structure

```text
covid-mpro-phytochemical-docking/

README.md
main.py
LICENSE
.gitignore

data/
    ligands.csv
    docking_results.csv

protocol/
    project_workflow.md

src/
    docking_analysis.py
```

## Selected Target

- Target protein: SARS-CoV-2 Main Protease
- PDB ID: 6LU7
- Database: RCSB Protein Data Bank

## Selected Ligands

- Quercetin
- Curcumin
- Kaempferol
- Andrographolide
- Berberine

## Workflow Summary

```text
Literature Review
        ↓
Target Protein Selection
        ↓
Ligand Selection
        ↓
Docking Result Organization
        ↓
Binding Affinity Ranking
        ↓
Summary Generation
```

## Technologies and Tools

- Python 3
- CSV data handling
- Object-Oriented Programming
- AutoDock Vina
- AutoDockTools
- PubChem
- RCSB Protein Data Bank
- PyMOL / Discovery Studio Visualizer

## Scientific Note

Molecular docking is a computational screening method used to estimate possible interactions between ligands and target proteins. Docking scores are preliminary and should not be interpreted as proof of biological activity.

Any compound identified through docking requires further validation through experimental, pharmacological, and clinical studies.

## Future Development

Planned extensions include:

- Integration of real docking log files
- Protein-ligand interaction tables
- Visualization of binding affinity rankings
- Addition of docking pose files
- PyMOL interaction snapshots
- Automated report generation
- ADMET screening workflow

## Author

**Ishitta Sarkar**

B.Tech Biotechnology

Areas of Interest:

- Molecular Docking
- Computational Drug Discovery
- Bioinformatics
- Computational Biology
- Precision Medicine
- Biomedical Data Science
