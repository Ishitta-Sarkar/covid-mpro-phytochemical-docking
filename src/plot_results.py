"""Write ``results/vina_scores.png`` from the results table."""

from mpro_docking.paths import data_dir, results_dir
from mpro_docking.plotting import plot_vina_scores


def main() -> None:
    plot_vina_scores(
        data_dir() / "docking_results.csv",
        results_dir() / "vina_scores.png",
        title="SARS-CoV-2 Mpro (PDB 6LU7)",
    )
    other = data_dir() / "docking_results_6w63.csv"
    if other.exists():
        plot_vina_scores(
            other,
            results_dir() / "vina_scores_6w63.png",
            title="SARS-CoV-2 Mpro (PDB 6W63)",
        )


if __name__ == "__main__":
    main()
