"""Write ``results/vina_scores.png`` from the results table."""

from mpro_docking.paths import data_dir, results_dir
from mpro_docking.plotting import plot_vina_scores


def main() -> None:
    plot_vina_scores(data_dir() / "docking_results.csv", results_dir() / "vina_scores.png")


if __name__ == "__main__":
    main()
