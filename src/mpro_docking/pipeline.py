"""End-to-end fetch, preparation, docking, redocking and tables."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import platform
import urllib.error
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from mpro_docking import __version__
from mpro_docking.chemistry import assemble_crystal_ligand, residue_atom_key
from mpro_docking.constants import (
    CPU,
    EXHAUSTIVENESS,
    PDB_ID,
    RECEPTOR_CHAIN,
    SEED,
)
from mpro_docking.fetch import fetch_pdb, fetch_pubchem_sdf, load_ligand_sdf, sha256_file
from mpro_docking.geometry import cubic_search_box, distance
from mpro_docking.interactions import (
    contact_residues,
    hydrogen_bonds,
    min_intermolecular_distance,
    parse_pdbqt_atoms,
    write_contacts,
    write_hydrogen_bonds,
)
from mpro_docking.paths import repo_root
from mpro_docking.pdbio import parse_pdb
from mpro_docking.plotting import plot_vina_scores
from mpro_docking.prepare import (
    obabel_pdb,
    prepare_ligand_pdbqt,
    prepare_receptor_pdbqt,
    write_receptor_pdb,
    write_sdf,
)
from mpro_docking.rmsd import (
    heavy_atom_rmsd,
    mol_with_coordinates,
    parse_index_map,
    parse_pdbqt_coordinates,
)
from mpro_docking.vina_io import (
    format_vina_config,
    parse_score_only,
    parse_vina_modes,
    run_vina,
    torsdof,
)

LOGGER = logging.getLogger("mpro_docking")

RESULT_FIELDS = [
    "rank",
    "ligand_name",
    "pubchem_cid",
    "compound_class",
    "pdb_id",
    "vina_score_kcal_mol",
    "torsdof",
    "seed",
    "exhaustiveness",
    "cpu",
    "log_file",
    "pose_file",
]


def _read_ligands(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 5:
        raise ValueError(f"{path} should list five ligands, found {len(rows)}")
    return rows


def _relative(path: Path, root: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def _dock_one(
    root: Path,
    config: Path,
    ligand_pdbqt: Path,
    pose: Path,
    log: Path,
) -> tuple[float, int]:
    run_vina(config, ligand_pdbqt, pose, log, cpu=CPU, cwd=root)
    pose_modes = parse_vina_modes(pose.read_text(encoding="utf-8"))
    log_modes = parse_vina_modes(log.read_text(encoding="utf-8"))
    if not pose_modes:
        raise RuntimeError(f"no Vina scores in {pose}")
    if log_modes and abs(log_modes[0].score - pose_modes[0].score) > 0.011:
        raise RuntimeError(
            f"log and pose scores disagree for {ligand_pdbqt.name}: "
            f"{log_modes[0].score} vs {pose_modes[0].score}"
        )
    log_text = log.read_text(encoding="utf-8")
    if f"random seed: {SEED}" not in log_text and f"seed {SEED}" not in log_text.lower():
        LOGGER.warning("Vina log does not echo seed %s (%s)", SEED, log.name)
    return pose_modes[0].score, torsdof(ligand_pdbqt.read_text(encoding="utf-8"))


def _roundtrip_rmsd(reference, pdbqt: Path) -> float:
    text = pdbqt.read_text(encoding="utf-8")
    mapped = mol_with_coordinates(
        reference,
        parse_pdbqt_coordinates(text),
        parse_index_map(text),
    )
    return heavy_atom_rmsd(mapped, reference).direct


def _pose_mol(reference, input_pdbqt: Path, output_pdbqt: Path, model: int = 1):
    mapped = mol_with_coordinates(
        reference,
        parse_pdbqt_coordinates(output_pdbqt.read_text(encoding="utf-8"), model=model),
        parse_index_map(input_pdbqt.read_text(encoding="utf-8")),
    )
    return mapped


def _atom_xyz(mol, resname: str, resseq: int, name: str) -> tuple[float, float, float]:
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        key = residue_atom_key(atom)
        if key == (resname, resseq, name):
            pos = conf.GetAtomPosition(atom.GetIdx())
            return (pos.x, pos.y, pos.z)
    raise ValueError(f"atom {resname} {resseq} {name} not found")


def run(refresh: bool = False) -> dict:
    root = repo_root()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    raw = root / "structures" / "raw"
    receptor_dir = root / "structures" / "receptor"
    ligand_dir = root / "structures" / "ligands"
    pose_dir = root / "structures" / "poses"
    log_dir = root / "docking" / "logs"
    result_dir = root / "results"

    LOGGER.info("Fetching PDB %s", PDB_ID)
    try:
        pdb_path = fetch_pdb(PDB_ID, raw / f"{PDB_ID}.pdb", refresh=refresh)
    except urllib.error.URLError as exc:
        raise RuntimeError(f"could not download PDB {PDB_ID}: {exc}") from exc

    LOGGER.info("Preparing receptor chain %s", RECEPTOR_CHAIN)
    clean_pdb = receptor_dir / f"{PDB_ID}_chainA.pdb"
    n_protein = write_receptor_pdb(pdb_path, clean_pdb)
    receptor_pdbqt = receptor_dir / f"{PDB_ID}_receptor.pdbqt"
    prepare_receptor_pdbqt(clean_pdb, receptor_pdbqt)

    LOGGER.info("Building co-crystallised inhibitor from chain C")
    crystal = assemble_crystal_ligand(pdb_path, root / "data" / "ccd")
    crystal_sdf = ligand_dir / "n3_crystal.sdf"
    write_sdf(crystal, crystal_sdf)
    obabel_pdb(crystal_sdf, ligand_dir / "n3_crystal.pdb")
    n3_pdbqt = ligand_dir / "n3.pdbqt"
    prepare_ligand_pdbqt(crystal, n3_pdbqt, "N3")
    roundtrip = _roundtrip_rmsd(crystal, n3_pdbqt)
    if roundtrip > 0.01:
        raise RuntimeError(
            f"N3 PDBQT does not preserve crystal heavy atoms (RMSD {roundtrip:.4f} A)"
        )

    heavy_coords = [
        crystal.GetConformer().GetAtomPosition(atom.GetIdx())
        for atom in crystal.GetAtoms()
        if atom.GetAtomicNum() != 1
    ]
    heavy_xyz = [(p.x, p.y, p.z) for p in heavy_coords]
    center, size = cubic_search_box(heavy_xyz)
    config_path = root / "docking" / "vina_config.txt"
    config_path.write_text(
        format_vina_config(_relative(receptor_pdbqt, root), center, size),
        encoding="utf-8",
    )
    LOGGER.info(
        "Search box centre (%.3f, %.3f, %.3f), size %.0f A",
        center[0],
        center[1],
        center[2],
        size[0],
    )

    LOGGER.info("Scoring the crystallographic N3 pose without a search")
    score_log = log_dir / "n3_score_only.log"
    score_text = run_vina(
        config_path,
        n3_pdbqt,
        None,
        score_log,
        extra_args=["--score_only"],
        cwd=root,
    )
    crystal_score = parse_score_only(score_text)

    LOGGER.info("Redocking N3 (exhaustiveness %s, seed %s, cpu %s)", EXHAUSTIVENESS, SEED, CPU)
    n3_pose = pose_dir / "n3_out.pdbqt"
    n3_score, n3_tors = _dock_one(root, config_path, n3_pdbqt, n3_pose, log_dir / "n3.log")
    docked = _pose_mol(crystal, n3_pdbqt, n3_pose)
    rmsd = heavy_atom_rmsd(docked, crystal)
    docked_sdf = pose_dir / "n3_top.sdf"
    write_sdf(docked, docked_sdf)
    obabel_pdb(docked_sdf, pose_dir / "n3_top.pdb")
    LOGGER.info(
        "Redocking heavy-atom RMSD (symmetry-corrected, no superposition): %.3f A",
        rmsd.symmetry,
    )
    LOGGER.info("Direct heavy-atom RMSD: %.3f A", rmsd.direct)

    receptor_atoms = parse_pdbqt_atoms(receptor_pdbqt)
    sg = next(
        atom.xyz
        for atom in receptor_atoms
        if atom.resname == "CYS" and atom.resseq == 145 and atom.name == "SG"
    )
    crystal_c20 = _atom_xyz(crystal, "PJE", 5, "C20")
    docked_c20 = _atom_xyz(docked, "PJE", 5, "C20")
    n3_pose_atoms = parse_pdbqt_atoms(n3_pose)

    ligand_rows = _read_ligands(root / "data" / "ligands.csv")
    provenance = []
    result_rows = []
    all_hbonds = []
    all_contacts = []

    for row in ligand_rows:
        name = row["ligand_name"]
        cid = row["pubchem_cid"]
        slug = name.lower()
        LOGGER.info("Fetching and docking %s (PubChem CID %s)", name, cid)
        sdf_path, source = fetch_pubchem_sdf(cid, raw / f"{slug}_{cid}.sdf", refresh=refresh)
        mol, source = load_ligand_sdf(sdf_path, source)
        mol.SetProp("_Name", name)
        prepared_sdf = ligand_dir / f"{slug}.sdf"
        write_sdf(mol, prepared_sdf)
        obabel_pdb(prepared_sdf, ligand_dir / f"{slug}.pdb")
        pdbqt = ligand_dir / f"{slug}.pdbqt"
        prepare_ligand_pdbqt(mol, pdbqt, name)
        trip = _roundtrip_rmsd(mol, pdbqt)
        if trip > 0.05:
            raise RuntimeError(f"{name} PDBQT moved heavy atoms by {trip:.3f} A before docking")
        pose = pose_dir / f"{slug}_out.pdbqt"
        score, n_tors = _dock_one(root, config_path, pdbqt, pose, log_dir / f"{slug}.log")
        pose_mol = _pose_mol(mol, pdbqt, pose)
        write_sdf(pose_mol, pose_dir / f"{slug}_top.sdf")
        ligand_atoms = parse_pdbqt_atoms(pose)
        all_hbonds.extend(hydrogen_bonds(name, receptor_atoms, ligand_atoms))
        all_contacts.extend(contact_residues(name, receptor_atoms, ligand_atoms))
        result_rows.append(
            {
                "ligand_name": name,
                "pubchem_cid": cid,
                "compound_class": row["compound_class"],
                "pdb_id": PDB_ID,
                "vina_score_kcal_mol": f"{score:.3f}",
                "torsdof": str(n_tors),
                "seed": str(SEED),
                "exhaustiveness": str(EXHAUSTIVENESS),
                "cpu": str(CPU),
                "log_file": _relative(log_dir / f"{slug}.log", root),
                "pose_file": _relative(pose, root),
            }
        )
        provenance.append(
            {
                "ligand_name": name,
                "pubchem_cid": cid,
                "source": source,
                "sdf_sha256": sha256_file(sdf_path),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "formal_charge": Chem.GetFormalCharge(mol),
                "heavy_atoms": mol.GetNumHeavyAtoms(),
            }
        )
        LOGGER.info("%s Vina score %.3f kcal/mol (%s torsions)", name, score, n_tors)

    result_rows.sort(key=lambda item: (float(item["vina_score_kcal_mol"]), item["ligand_name"]))
    for rank, row in enumerate(result_rows, start=1):
        row["rank"] = str(rank)

    results_csv = root / "data" / "docking_results.csv"
    with results_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(result_rows)

    top_csv = result_dir / "top_candidates.csv"
    with top_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "ligand_name", "vina_score_kcal_mol", "note"],
        )
        writer.writeheader()
        for row in result_rows:
            writer.writerow(
                {
                    "rank": row["rank"],
                    "ligand_name": row["ligand_name"],
                    "vina_score_kcal_mol": row["vina_score_kcal_mol"],
                    "note": "Ranked by Vina score in this run. Not a measured affinity.",
                }
            )

    protein_atoms, _links = parse_pdb(pdb_path)
    his41 = next(
        atom.xyz
        for atom in protein_atoms
        if atom.chain == "A" and atom.resname == "HIS" and atom.resseq == 41 and atom.name == "CA"
    )
    redock = {
        "ligand": "N3",
        "description": (
            "Chain C of PDB 6LU7 (02J-ALA-VAL-LEU-PJE-010), Cys145 SG-PJE C20 bond omitted"
        ),
        "pdb_id": PDB_ID,
        "heavy_atom_rmsd_symmetry_A": round(rmsd.symmetry, 3),
        "heavy_atom_rmsd_direct_A": round(rmsd.direct, 3),
        "n_heavy_atoms": rmsd.n_heavy,
        "n_automorphisms": rmsd.n_automorphisms,
        "vina_score_kcal_mol": round(n3_score, 3),
        "crystal_pose_score_kcal_mol": None if crystal_score is None else round(crystal_score, 3),
        "torsdof": n3_tors,
        "seed": SEED,
        "exhaustiveness": EXHAUSTIVENESS,
        "cpu": CPU,
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": size[0],
        "size_y": size[1],
        "size_z": size[2],
        "crystal_c20_sg_distance_A": round(distance(crystal_c20, sg), 3),
        "docked_c20_sg_distance_A": round(distance(docked_c20, sg), 3),
        "docked_min_heavy_distance_A": round(
            min_intermolecular_distance(receptor_atoms, n3_pose_atoms), 3
        ),
        "formula": rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(crystal)),
        "protein_atoms_chain_A": n_protein,
    }
    redock_path = result_dir / "redocking_n3.json"
    redock_path.write_text(json.dumps(redock, indent=2) + "\n", encoding="utf-8")
    with (result_dir / "redocking_n3.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(redock.keys()))
        writer.writeheader()
        writer.writerow(redock)

    write_hydrogen_bonds(all_hbonds, result_dir / "hydrogen_bonds.csv")
    write_contacts(all_contacts, result_dir / "residue_contacts.csv")
    figure = result_dir / "vina_scores.png"
    plot_vina_scores(results_csv, figure)

    metadata = {
        "software_version": __version__,
        "date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "vina_package": "1.2.7",
        "vina_binary": _vina_version(),
        "meeko": version("meeko"),
        "rdkit": Chem.rdBase.rdkitVersion,
        "openbabel": _command_version(["obabel", "-V"]),
        "seed": SEED,
        "exhaustiveness": EXHAUSTIVENESS,
        "cpu": CPU,
        "box_center": list(center),
        "box_size": list(size),
        "box_rule": (
            "Cubic edge is the smallest even angstrom length that is at least 22 "
            "and at least twice the sum of 4 A padding and the maximum absolute "
            "coordinate deviation of the 49 crystallographic heavy atoms from their centroid. "
            "The centre written to the config is that centroid rounded to 3 decimals."
        ),
        "his41_ca": [round(v, 3) for v in his41],
        "distance_center_to_his41_ca_A": round(distance(center, his41), 3),
        "cys145_sg": [round(v, 3) for v in sg],
        "distance_center_to_cys145_sg_A": round(distance(center, sg), 3),
        "pdb_sha256": sha256_file(pdb_path),
        "receptor_pdbqt": _relative(receptor_pdbqt, root),
        "ligands": provenance,
        "redocking_rmsd_symmetry_A": redock["heavy_atom_rmsd_symmetry_A"],
        "scores": [
            {
                "ligand_name": row["ligand_name"],
                "vina_score_kcal_mol": float(row["vina_score_kcal_mol"]),
            }
            for row in result_rows
        ],
    }
    (result_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    _write_pymol(root)
    LOGGER.info("Wrote %s", _relative(results_csv, root))
    LOGGER.info("Wrote %s", _relative(figure, root))
    return {"redocking": redock, "results": result_rows, "metadata": metadata}


def _vina_version() -> str:
    return _command_version(["vina", "--version"])


def _command_version(cmd: list[str]) -> str:
    import subprocess

    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
    text = (completed.stdout or completed.stderr).strip().splitlines()
    return text[0] if text else "unknown"


def _write_pymol(root: Path) -> None:
    script = root / "docking" / "overlay_n3.pml"
    script.write_text(
        "\n".join(
            [
                "reinitialize",
                "bg_color white",
                "set ray_opaque_background, 1",
                "load structures/receptor/6LU7_chainA.pdb, receptor",
                "load structures/ligands/n3_crystal.pdb, crystal",
                "load structures/poses/n3_top.pdb, redocked",
                "hide everything",
                "hide (elem H)",
                "show cartoon, receptor",
                "show sticks, crystal or redocked",
                "set stick_radius, 0.18",
                "color grey70, receptor",
                "color cyan, crystal",
                "color salmon, redocked",
                "select site, chain A and (resi 41 or resi 145)",
                "show sticks, site",
                "color yellow, site",
                "orient crystal",
                "zoom crystal, 8",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Dock five phytochemicals against SARS-CoV-2 Mpro (PDB 6LU7) with AutoDock Vina."
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Download the PDB and PubChem SDFs again even if local copies exist.",
    )
    args = parser.parse_args(argv)
    run(refresh=args.refresh)
    return 0
