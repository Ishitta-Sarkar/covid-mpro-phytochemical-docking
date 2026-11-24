"""Hydrogen bonds and first-shell contacts from prepared PDBQT coordinates.

Hydrogens are those already present in the Meeko receptor and the docked
ligand. A hydrogen bond is accepted when the H...A distance is at most
2.5 Å and the D-H...A angle is at least 120 degrees (Baker and Hubbard).
Donors and acceptors are N and O. A contact residue has any heavy atom
within 4.0 Å of a ligand heavy atom. These are geometric filters, not a
substitute for an interaction-energy calculation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from mpro_docking.constants import CONTACT_CUTOFF, HBOND_MAX_HA, HBOND_MIN_ANGLE
from mpro_docking.geometry import distance, is_hydrogen_bond

_ACCEPTOR_TYPES = {"OA", "NA", "N", "OS", "OD", "ND"}
_POLAR_H = {"HD"}
_HEAVY_ELEMENTS = {"C", "N", "O", "S", "F", "P", "CL", "BR", "I"}


@dataclass(frozen=True)
class PdbqtAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    xyz: tuple[float, float, float]
    atom_type: str

    @property
    def element(self) -> str:
        if self.atom_type in {"HD", "H"}:
            return "H"
        if self.atom_type in {"OA", "OS", "OD"}:
            return "O"
        if self.atom_type in {"NA", "N", "ND"}:
            return "N"
        if self.atom_type in {"SA", "S"}:
            return "S"
        if self.atom_type in {"CL", "BR"}:
            return self.atom_type.title()
        if self.atom_type in {"A", "C"}:
            return "C"
        if self.atom_type in {"F", "P", "I"}:
            return self.atom_type
        letters = "".join(ch for ch in self.name if ch.isalpha())
        return letters[:1].upper()


def parse_pdbqt_atoms(path: Path, model: int = 1) -> list[PdbqtAtom]:
    models: list[list[PdbqtAtom]] = []
    current: list[PdbqtAtom] | None = None
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            current = []
            models.append(current)
            continue
        if line.startswith("ENDMDL"):
            current = None
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if saw_model and current is None:
            continue
        if current is None:
            current = []
            models.append(current)
        tokens = line.split()
        current.append(
            PdbqtAtom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resname=line[17:20].strip() or "UNL",
                chain=line[21:22].strip() or "L",
                resseq=int((line[22:26].strip() or "0")),
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                atom_type=tokens[-1],
            )
        )
    if not models:
        raise ValueError(f"no atoms in {path}")
    return models[model - 1]


@dataclass(frozen=True)
class HydrogenBond:
    ligand_name: str
    donor_side: str
    donor_label: str
    acceptor_label: str
    distance_a: float
    angle_deg: float


def hydrogen_bonds(
    ligand_name: str,
    receptor_atoms: list[PdbqtAtom],
    ligand_atoms: list[PdbqtAtom],
) -> list[HydrogenBond]:
    found: list[HydrogenBond] = []
    found.extend(_pairs(ligand_name, "ligand", ligand_atoms, receptor_atoms))
    found.extend(_pairs(ligand_name, "protein", receptor_atoms, ligand_atoms))
    found.sort(key=lambda row: (row.distance_a, row.donor_label, row.acceptor_label))
    return found


def _pairs(
    ligand_name: str,
    donor_side: str,
    donors: list[PdbqtAtom],
    acceptors: list[PdbqtAtom],
) -> list[HydrogenBond]:
    heavies = [atom for atom in donors if atom.element in _HEAVY_ELEMENTS]
    acceptor_atoms = [
        atom
        for atom in acceptors
        if atom.element in {"N", "O"} or atom.atom_type in _ACCEPTOR_TYPES
    ]
    bonds: dict[tuple[str, str], HydrogenBond] = {}
    for hydrogen in donors:
        if hydrogen.atom_type not in _POLAR_H and hydrogen.element != "H":
            continue
        if hydrogen.atom_type not in _POLAR_H:
            continue
        parent = _parent(hydrogen, heavies)
        if parent is None or parent.element not in {"N", "O"}:
            continue
        for acceptor in acceptor_atoms:
            ok, dist, ang = is_hydrogen_bond(
                parent.xyz,
                hydrogen.xyz,
                acceptor.xyz,
                max_ha=HBOND_MAX_HA,
                min_angle=HBOND_MIN_ANGLE,
            )
            if not ok:
                continue
            donor_label = _label(parent)
            acceptor_label = _label(acceptor)
            key = (donor_label, acceptor_label)
            row = HydrogenBond(ligand_name, donor_side, donor_label, acceptor_label, dist, ang)
            previous = bonds.get(key)
            if previous is None or dist < previous.distance_a:
                bonds[key] = row
    return list(bonds.values())


def _parent(hydrogen: PdbqtAtom, heavies: list[PdbqtAtom]) -> PdbqtAtom | None:
    best = None
    best_d = 1.25
    for atom in heavies:
        dist = distance(hydrogen.xyz, atom.xyz)
        if dist < best_d:
            best = atom
            best_d = dist
    return best


def _label(atom: PdbqtAtom) -> str:
    return f"{atom.chain}:{atom.resname}{atom.resseq}:{atom.name}"


def contact_residues(
    ligand_name: str,
    receptor_atoms: list[PdbqtAtom],
    ligand_atoms: list[PdbqtAtom],
    cutoff: float = CONTACT_CUTOFF,
) -> list[dict[str, str | int | float]]:
    ligand_heavy = [atom for atom in ligand_atoms if atom.element != "H"]
    closest: dict[tuple[str, str, int], float] = {}
    for protein in receptor_atoms:
        if protein.element == "H":
            continue
        for ligand in ligand_heavy:
            dist = distance(protein.xyz, ligand.xyz)
            if dist > cutoff:
                continue
            key = (protein.chain, protein.resname, protein.resseq)
            previous = closest.get(key)
            if previous is None or dist < previous:
                closest[key] = dist
    rows = [
        {
            "ligand_name": ligand_name,
            "chain": chain,
            "resname": resname,
            "resseq": resseq,
            "min_distance_A": round(dist, 3),
        }
        for (chain, resname, resseq), dist in closest.items()
    ]
    rows.sort(
        key=lambda row: (float(row["min_distance_A"]), str(row["resname"]), int(row["resseq"]))
    )
    return rows


def write_hydrogen_bonds(rows: list[HydrogenBond], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ligand_name",
                "donor_side",
                "donor",
                "acceptor",
                "h_a_distance_A",
                "d_h_a_angle_deg",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "ligand_name": row.ligand_name,
                    "donor_side": row.donor_side,
                    "donor": row.donor_label,
                    "acceptor": row.acceptor_label,
                    "h_a_distance_A": f"{row.distance_a:.3f}",
                    "d_h_a_angle_deg": f"{row.angle_deg:.1f}",
                }
            )


def write_contacts(rows: list[dict[str, str | int | float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ligand_name", "chain", "resname", "resseq", "min_distance_A"],
        )
        writer.writeheader()
        writer.writerows(rows)


def min_intermolecular_distance(
    receptor_atoms: list[PdbqtAtom], ligand_atoms: list[PdbqtAtom]
) -> float:
    best = None
    for protein in receptor_atoms:
        if protein.element == "H":
            continue
        for ligand in ligand_atoms:
            if ligand.element == "H":
                continue
            dist = distance(protein.xyz, ligand.xyz)
            if best is None or dist < best:
                best = dist
    if best is None:
        raise ValueError("no heavy atoms to compare")
    return best
