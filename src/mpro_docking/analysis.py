"""Rank a Vina results table. Lower scores come first."""

from __future__ import annotations

import csv
from pathlib import Path

SCORE_COLUMN = "vina_score_kcal_mol"


class DockingAnalysis:
    def __init__(self, ligand_file: str | Path, result_file: str | Path):
        self.ligands = self.load_csv(ligand_file)
        self.results = self.load_csv(result_file)

    @staticmethod
    def load_csv(file_path: str | Path) -> list[dict[str, str]]:
        with Path(file_path).open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def rank_ligands(self) -> list[dict[str, str]]:
        return sorted(
            self.results,
            key=lambda row: (float(row[SCORE_COLUMN]), row["ligand_name"]),
        )

    def best_candidate(self) -> dict[str, str]:
        ranked = self.rank_ligands()
        if not ranked:
            raise ValueError("no docking results to rank")
        return ranked[0]

    def generate_summary(self) -> str:
        lines = [
            "Virtual screening summary",
            "Scores are AutoDock Vina docking scores (kcal/mol).",
            "They are not measured binding affinities.",
            "",
        ]
        for rank, row in enumerate(self.rank_ligands(), start=1):
            lines.append(f"{rank}. {row['ligand_name']} | {row[SCORE_COLUMN]} kcal/mol")
        return "\n".join(lines) + "\n"
