"""Prepare a receptor PDBQT and ligand PDBQT files for Vina."""

from __future__ import annotations

import subprocess
from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

from mpro_docking.constants import RECEPTOR_CHAIN
from mpro_docking.pdbio import parse_pdb, write_atom_lines
from mpro_docking.rmsd import parse_index_map


def write_receptor_pdb(source_pdb: Path, dest_pdb: Path, chain: str = RECEPTOR_CHAIN) -> int:
    """Write protein ATOM records for one chain, using the deposited lines."""

    atoms, _links = parse_pdb(source_pdb)
    lines = [
        atom.line
        for atom in atoms
        if atom.record == "ATOM" and atom.chain == chain and atom.element != "H"
    ]
    if not lines:
        raise ValueError(f"no protein atoms for chain {chain} in {source_pdb}")
    write_atom_lines(lines, dest_pdb)
    return len(lines)


def prepare_receptor_pdbqt(clean_pdb: Path, pdbqt: Path) -> None:
    """Add template hydrogens and Gasteiger charges with Meeko and write a PDBQT."""

    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    basename = pdbqt.with_suffix("")
    cmd = [
        "mk_prepare_receptor.py",
        "--read_pdb",
        str(clean_pdb),
        "-o",
        str(basename),
        "-p",
        str(pdbqt),
        "--default_altloc",
        "A",
    ]
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0 or not pdbqt.exists():
        raise RuntimeError(
            "mk_prepare_receptor.py failed:\n" + (completed.stderr or completed.stdout)[-2000:]
        )


def prepare_ligand_pdbqt(mol: Chem.Mol, pdbqt: Path, name: str) -> dict[int, int]:
    """Write a Meeko PDBQT and return the RDKit-index to serial map."""

    named = Chem.Mol(mol)
    named.SetProp("_Name", name)
    preparer = MoleculePreparation()
    setups = preparer.prepare(named)
    if len(setups) != 1:
        raise RuntimeError(f"{name}: Meeko returned {len(setups)} setups; expected one")
    pdbqt_text, ok, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not ok:
        pdbqt_text, ok, error = PDBQTWriterLegacy.write_string(
            setups[0], add_index_map=True, bad_charge_ok=True
        )
    if not ok:
        raise RuntimeError(f"{name}: Meeko could not write a PDBQT: {error}")
    pdbqt.parent.mkdir(parents=True, exist_ok=True)
    pdbqt.write_text(pdbqt_text, encoding="utf-8")
    return parse_index_map(pdbqt_text)


def write_sdf(mol: Chem.Mol, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    try:
        writer.write(mol)
    finally:
        writer.close()


def obabel_pdb(sdf: Path, pdb: Path) -> None:
    """Write a PDB copy with Open Babel so the ligand can be opened outside RDKit."""

    pdb.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        ["obabel", str(sdf), "-O", str(pdb)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"obabel failed for {sdf.name}: {completed.stderr[-1000:]}")
