"""Print the ranked Vina scores from the repository root."""

from pathlib import Path

from docking_analysis import DockingAnalysis

ROOT = Path(__file__).resolve().parents[1]


def _print_table(title: str, path) -> None:
    analysis = DockingAnalysis(ROOT / "data" / "ligands.csv", path)
    print(title)
    print("Vina scores are not measured binding affinities.")
    for rank, row in enumerate(analysis.rank_ligands(), start=1):
        print(f"{rank}. {row['ligand_name']} ({row['vina_score_kcal_mol']} kcal/mol)")
    print()


def main() -> None:
    _print_table(
        "PDB 6LU7 (covalent complex; not the pose-reproduction control)",
        ROOT / "data" / "docking_results.csv",
    )
    other = ROOT / "data" / "docking_results_6w63.csv"
    if other.exists():
        _print_table(
            "PDB 6W63 (non-covalent complex used to validate the pose protocol)",
            other,
        )


if __name__ == "__main__":
    main()
