import csv
from pathlib import Path

from mpro_docking.analysis import DockingAnalysis
from mpro_docking.paths import data_dir
from mpro_docking.plotting import plot_vina_scores


def test_ranking_uses_more_negative_scores_first(tmp_path: Path):
    ligands = tmp_path / "ligands.csv"
    results = tmp_path / "results.csv"
    ligands.write_text(
        "ligand_name,pubchem_cid,compound_class\nA,1,x\nB,2,y\n",
        encoding="utf-8",
    )
    results.write_text(
        "ligand_name,vina_score_kcal_mol\nA,-6.0\nB,-8.5\n",
        encoding="utf-8",
    )
    analysis = DockingAnalysis(ligands, results)
    ranked = analysis.rank_ligands()
    assert [row["ligand_name"] for row in ranked] == ["B", "A"]
    assert analysis.best_candidate()["ligand_name"] == "B"
    summary = analysis.generate_summary()
    assert "not measured binding affinities" in summary
    assert "1. B |" in summary
    assert "2. A |" in summary


def test_ligand_table_has_the_five_pubchem_records():
    with (data_dir() / "ligands.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [(row["ligand_name"], row["pubchem_cid"]) for row in rows] == [
        ("Quercetin", "5280343"),
        ("Curcumin", "969516"),
        ("Kaempferol", "5280863"),
        ("Andrographolide", "5318517"),
        ("Berberine", "2353"),
    ]


def test_plot_writes_png(tmp_path: Path):
    results = tmp_path / "results.csv"
    results.write_text(
        "ligand_name,vina_score_kcal_mol\nQuercetin,-7.100\nCurcumin,-6.200\n",
        encoding="utf-8",
    )
    figure = tmp_path / "vina_scores.png"
    plot_vina_scores(results, figure)
    assert figure.read_bytes().startswith(b"\x89PNG")
