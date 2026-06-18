from src.docking_analysis import DockingAnalysis


analysis = DockingAnalysis(
    "data/ligands.csv",
    "data/docking_results.csv"
)

print("COVID-19 Mpro Docking Analysis")
print("=" * 40)

print("\nRanked Ligands")
print("-" * 40)

for rank, ligand in enumerate(analysis.rank_ligands(), start=1):
    print(
        f"{rank}. {ligand['ligand_name']} "
        f"({ligand['binding_affinity_kcal_mol']} kcal/mol)"
    )

print("\nTop Candidate")
print("-" * 40)

best = analysis.best_candidate()

print("Ligand:", best["ligand_name"])
print("Binding Affinity:", best["binding_affinity_kcal_mol"], "kcal/mol")
print("Target:", best["target_protein"])

print("\nSummary")
print("-" * 40)

print(analysis.generate_summary())
