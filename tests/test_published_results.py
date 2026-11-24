"""Checks against the files produced by the Vina run committed in this repository."""

import csv
import json
from pathlib import Path

from mpro_docking.constants import EXHAUSTIVENESS, SEED
from mpro_docking.vina_io import parse_vina_modes

ROOT = Path(__file__).resolve().parents[1]


def _results():
    with (ROOT / "data" / "docking_results.csv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_results_table_matches_the_ligand_list_and_the_logs():
    with (ROOT / "data" / "ligands.csv").open(encoding="utf-8", newline="") as handle:
        ligands = list(csv.DictReader(handle))
    rows = _results()
    assert len(rows) == 5
    assert {row["ligand_name"] for row in rows} == {row["ligand_name"] for row in ligands}
    assert [(row["ligand_name"], row["vina_score_kcal_mol"]) for row in rows] == [
        ("Kaempferol", "-7.745"),
        ("Berberine", "-7.679"),
        ("Quercetin", "-7.518"),
        ("Curcumin", "-6.910"),
        ("Andrographolide", "-6.883"),
    ]
    for row in rows:
        score = float(row["vina_score_kcal_mol"])
        assert score < 0
        assert int(row["seed"]) == SEED
        assert int(row["exhaustiveness"]) == EXHAUSTIVENESS
        log_modes = parse_vina_modes((ROOT / row["log_file"]).read_text(encoding="utf-8"))
        pose_modes = parse_vina_modes((ROOT / row["pose_file"]).read_text(encoding="utf-8"))
        assert pose_modes[0].score == pytest_approx(score)
        assert log_modes[0].score == pytest_approx(score)


def test_redocking_record_and_figure():
    redock = json.loads((ROOT / "results" / "redocking_n3.json").read_text(encoding="utf-8"))
    assert redock["heavy_atom_rmsd_symmetry_A"] == 7.266
    assert redock["heavy_atom_rmsd_direct_A"] == 7.295
    assert redock["n_heavy_atoms"] == 49
    assert redock["crystal_pose_score_kcal_mol"] == -3.47
    assert redock["vina_score_kcal_mol"] == -8.52
    assert (ROOT / "results" / "vina_scores.png").read_bytes().startswith(b"\x89PNG")
    assert (ROOT / "results" / "n3_redock_overlay.png").read_bytes().startswith(b"\x89PNG")
    config = (ROOT / "docking" / "vina_config.txt").read_text(encoding="utf-8")
    assert f"center_x = {redock['center_x']}" in config
    assert "seed = 42" in config
    assert "size_x = 30" in config


def pytest_approx(value: float):
    import pytest

    return pytest.approx(value, abs=0.001)
