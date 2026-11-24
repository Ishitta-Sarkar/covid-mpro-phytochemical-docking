"""Heavy-atom RMSD in a shared coordinate frame.

Docking RMSD must not superimpose the two poses. Symmetry correction uses
graph automorphisms so constitutionally identical atoms (phenyl flips,
gem-dimethyl groups) can be matched. The identity mapping is always included,
so the symmetry-corrected value cannot exceed the direct RMSD.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Geometry import Point3D


@dataclass(frozen=True)
class RmsdResult:
    symmetry: float
    direct: float
    n_heavy: int
    n_automorphisms: int


def parse_index_map(pdbqt_text: str) -> dict[int, int]:
    """Map RDKit atom index to PDBQT serial from Meeko ``REMARK INDEX MAP`` lines."""

    mapping: dict[int, int] = {}
    for line in pdbqt_text.splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue
        nums = [int(token) for token in line.split()[3:]]
        if len(nums) % 2:
            raise ValueError(f"odd number of INDEX MAP values: {line}")
        for rdkit_1based, serial in zip(nums[0::2], nums[1::2]):
            mapping[rdkit_1based - 1] = serial
    if not mapping:
        raise ValueError("PDBQT file has no REMARK INDEX MAP lines")
    return mapping


def parse_pdbqt_coordinates(
    pdbqt_text: str, model: int = 1
) -> dict[int, tuple[float, float, float]]:
    """Return serial -> xyz for one pose.

    Vina output is split into MODEL blocks. An input PDBQT has no MODEL
    records and is treated as a single pose.
    """

    models: list[dict[int, tuple[float, float, float]]] = []
    current: dict[int, tuple[float, float, float]] | None = None
    saw_model = False
    for line in pdbqt_text.splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            current = {}
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
            current = {}
            models.append(current)
        serial = int(line[6:11])
        current[serial] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    if not models:
        raise ValueError("PDBQT file has no atoms")
    if model < 1 or model > len(models):
        raise ValueError(f"model {model} is outside 1..{len(models)}")
    return models[model - 1]


def mol_with_coordinates(
    reference: Chem.Mol,
    serial_xyz: dict[int, tuple[float, float, float]],
    index_map: dict[int, int],
) -> Chem.Mol:
    """Copy ``reference`` and replace coordinates using a Meeko serial map."""

    out = Chem.Mol(reference)
    conf = out.GetConformer()
    n_atoms = reference.GetNumAtoms()
    missing = []
    for rdkit_index, serial in index_map.items():
        if rdkit_index < 0 or rdkit_index >= n_atoms:
            raise ValueError(
                f"index map points at RDKit atom {rdkit_index}, outside 0..{n_atoms - 1}"
            )
        if serial not in serial_xyz:
            missing.append(serial)
            continue
        x, y, z = serial_xyz[serial]
        conf.SetAtomPosition(rdkit_index, Point3D(x, y, z))
    if missing:
        raise ValueError(f"PDBQT is missing serials required by the index map: {missing[:8]}")
    # PDBQT is a united-atom format: nonpolar hydrogens are omitted.
    heavy_missing = [
        i
        for i, atom in enumerate(reference.GetAtoms())
        if atom.GetAtomicNum() != 1 and i not in index_map
    ]
    if heavy_missing:
        raise ValueError(f"index map misses heavy atoms: {heavy_missing[:8]}")
    return out


def heavy_atom_rmsd(probe: Chem.Mol, reference: Chem.Mol) -> RmsdResult:
    """RMSD of heavy atoms without superposition."""

    probe_h = Chem.RemoveHs(probe)
    ref_h = Chem.RemoveHs(reference)
    if probe_h.GetNumAtoms() != ref_h.GetNumAtoms():
        raise ValueError(
            f"heavy-atom counts differ: probe {probe_h.GetNumAtoms()}, "
            f"reference {ref_h.GetNumAtoms()}"
        )
    n = ref_h.GetNumAtoms()
    autos = list(
        ref_h.GetSubstructMatches(
            ref_h,
            useChirality=False,
            uniquify=False,
            maxMatches=10000,
        )
    )
    identity = tuple(range(n))
    if identity not in autos:
        autos.append(identity)
    probe_conf = probe_h.GetConformer()
    ref_conf = ref_h.GetConformer()
    direct = _mapped_rmsd(probe_conf, ref_conf, identity)
    best = min(_mapped_rmsd(probe_conf, ref_conf, mapping) for mapping in autos)
    return RmsdResult(symmetry=best, direct=direct, n_heavy=n, n_automorphisms=len(autos))


def _mapped_rmsd(probe_conf, ref_conf, mapping: tuple[int, ...]) -> float:
    total = 0.0
    for i, j in enumerate(mapping):
        p = probe_conf.GetAtomPosition(i)
        r = ref_conf.GetAtomPosition(j)
        total += (p.x - r.x) ** 2 + (p.y - r.y) ** 2 + (p.z - r.z) ** 2
    return math.sqrt(total / len(mapping))
