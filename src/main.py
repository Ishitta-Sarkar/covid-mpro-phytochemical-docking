"""Print the ranked Vina scores from the repository root."""

from pathlib import Path

from docking_analysis import DockingAnalysis

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    analysis = DockingAnalysis(
        ROOT / "data" / "ligands.csv",
        ROOT / "data" / "docking_results.csv",
    )
    print("SARS-CoV-2 Mpro docking scores")
    print("Vina scores are not measured binding affinities.")
    print()
    for rank, row in enumerate(analysis.rank_ligands(), start=1):
        print(f"{rank}. {row['ligand_name']} ({row['vina_score_kcal_mol']} kcal/mol)")
    print()
    best = analysis.best_candidate()
    print(f"Lowest score in this run: {best['ligand_name']}")
    print()
    print(analysis.generate_summary())


if __name__ == "__main__":
    main()
