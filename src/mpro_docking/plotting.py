"""Figure of the Vina scores from ``data/docking_results.csv``."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mpro_docking.analysis import SCORE_COLUMN, DockingAnalysis


def plot_vina_scores(results_csv: Path, figure_path: Path) -> None:
    ranked = sorted(
        DockingAnalysis(results_csv, results_csv).results,
        key=lambda row: (float(row[SCORE_COLUMN]), row["ligand_name"]),
    )
    names = [row["ligand_name"] for row in ranked]
    scores = [float(row[SCORE_COLUMN]) for row in ranked]
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    positions = list(range(len(names) - 1, -1, -1))
    ax.barh(positions, scores, color="#0072B2", height=0.62)
    ax.set_yticks(positions, names)
    ax.set_xlabel("Vina score (kcal/mol)")
    ax.set_title("SARS-CoV-2 Mpro (PDB 6LU7)")
    ax.axvline(0.0, color="#222222", linewidth=0.6)
    span = max(abs(score) for score in scores)
    ax.set_xlim(-span - 1.8, 0.4)
    for y, score in zip(positions, scores):
        ax.text(
            score - 0.12,
            y,
            f"{score:.3f}",
            va="center",
            ha="right",
            color="#222222",
            fontsize=9,
        )
    fig.text(
        0.125,
        0.02,
        "More negative scores rank higher in this run. A Vina score is not a measured binding affinity.",
        fontsize=8,
        color="#333333",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)
