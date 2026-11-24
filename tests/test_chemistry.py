import itertools
from collections import Counter
from pathlib import Path

import pytest
from rdkit import Chem

from mpro_docking.chemistry import (
    assemble_component_ligand,
    assemble_crystal_ligand,
    parse_ccd_bonds,
    residue_atom_key,
)
from mpro_docking.geometry import cubic_search_box, distance
from mpro_docking.pdbio import parse_pdb

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "6lu7_chainC.pdb"
CCD = ROOT / "data" / "ccd"


def _bond_pair(a, b):
    return tuple(sorted((a, b)))


EXPECTED_LINKS = {
    ("02J", 1, "C", "ALA", 2, "N"),
    ("LEU", 4, "C", "PJE", 5, "N"),
    ("PJE", 5, "C", "010", 6, "O"),
}
EXPECTED_INTERS = {
    _bond_pair(("02J", 1, "C"), ("ALA", 2, "N")),
    _bond_pair(("ALA", 2, "C"), ("VAL", 3, "N")),
    _bond_pair(("VAL", 3, "C"), ("LEU", 4, "N")),
    _bond_pair(("LEU", 4, "C"), ("PJE", 5, "N")),
    _bond_pair(("PJE", 5, "C"), ("010", 6, "O")),
}


def test_ccd_bond_orders():
    bonds = parse_ccd_bonds(CCD / "02J.cif")
    assert ("CA", "C4", "SING") in bonds
    assert ("C4", "C5", "DOUB") in bonds
    pje = parse_ccd_bonds(CCD / "PJE.cif")
    assert ("C21", "C20", "DOUB") in pje
    benzyl = parse_ccd_bonds(CCD / "010.cif")
    assert ("O", "C", "SING") in benzyl


def test_link_parser_keeps_chain_c_and_sees_covalent_record():
    _atoms, links = parse_pdb(FIXTURE)
    chain_c = {
        (link.resname1, link.resseq1, link.name1, link.resname2, link.resseq2, link.name2)
        for link in links
        if link.chain1 == "C" and link.chain2 == "C"
    }
    assert chain_c == EXPECTED_LINKS
    covalent = [link for link in links if link.resname1 == "CYS" and link.name1 == "SG"]
    assert len(covalent) == 1
    assert covalent[0].resname2 == "PJE"
    assert covalent[0].name2 == "C20"


def test_assembled_ligand_matches_crystal_heavy_atoms():
    mol = assemble_crystal_ligand(FIXTURE, CCD, relax_h=False)
    Chem.SanitizeMol(mol)
    heavy = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    assert len(heavy) == 49
    assert Chem.GetFormalCharge(mol) == 0
    assert not any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())

    pdb_atoms, _links = parse_pdb(FIXTURE)
    pdb_counts = Counter(atom.element for atom in pdb_atoms if atom.chain == "C")
    mol_counts = Counter(atom.GetSymbol() for atom in heavy)
    assert pdb_counts == mol_counts

    conf = mol.GetConformer()
    by_key = {}
    for atom in heavy:
        key = residue_atom_key(atom)
        pos = conf.GetAtomPosition(atom.GetIdx())
        by_key[key] = (pos.x, pos.y, pos.z)
    for atom in pdb_atoms:
        if atom.chain != "C":
            continue
        assert distance(by_key[atom.key], atom.xyz) < 1e-4

    warhead = None
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 1 or bond.GetEndAtom().GetAtomicNum() == 1:
            continue
        a = residue_atom_key(bond.GetBeginAtom())
        b = residue_atom_key(bond.GetEndAtom())
        if a[0] == "PJE" and b[0] == "PJE" and {a[2], b[2]} == {"C20", "C21"}:
            warhead = bond
    assert warhead is not None
    assert warhead.GetBondType() == Chem.BondType.DOUBLE

    inter = set()
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 1 or bond.GetEndAtom().GetAtomicNum() == 1:
            continue
        a = residue_atom_key(bond.GetBeginAtom())
        b = residue_atom_key(bond.GetEndAtom())
        if (a[0], a[1]) == (b[0], b[1]):
            continue
        inter.add(_bond_pair(a, b))
    assert inter == EXPECTED_INTERS

    coords = list(by_key.values())
    keys = list(by_key.keys())
    bonded = set()
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 1 or bond.GetEndAtom().GetAtomicNum() == 1:
            continue
        bonded.add(frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
    index = {residue_atom_key(atom): atom.GetIdx() for atom in heavy}
    close_missing = []
    for i, j in itertools.combinations(range(len(keys)), 2):
        dist = distance(coords[i], coords[j])
        pair = frozenset((index[keys[i]], index[keys[j]]))
        if dist < 1.65 and pair not in bonded:
            close_missing.append((keys[i], keys[j], round(dist, 3)))
        if pair in bonded and dist > 1.70:
            raise AssertionError(f"bonded pair is too long: {keys[i]} {keys[j]} {dist:.3f}")
    assert close_missing == []


def test_x77_is_a_complete_noncovalent_ligand():
    fixture = ROOT / "tests" / "fixtures" / "6w63_x77.pdb"
    mol = assemble_component_ligand(fixture, CCD / "X77.cif", "X77", "A", 401, "X77", relax_h=False)
    heavy = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    assert len(heavy) == 34
    assert Chem.GetFormalCharge(mol) == 0
    assert not any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    pdb_atoms, _links = parse_pdb(fixture)
    conf = mol.GetConformer()
    by_key = {}
    for atom in heavy:
        pos = conf.GetAtomPosition(atom.GetIdx())
        by_key[residue_atom_key(atom)] = (pos.x, pos.y, pos.z)
    for atom in pdb_atoms:
        assert distance(by_key[atom.key], atom.xyz) < 1e-4
    center, size = cubic_search_box([atom.xyz for atom in pdb_atoms])
    assert size == (22.0, 22.0, 22.0)
    assert center == (-20.096, 18.844, -27.394)


def test_covalent_link_is_rejected_for_the_noncovalent_builder():
    with pytest.raises(ValueError, match="covalent"):
        assemble_component_ligand(FIXTURE, CCD / "PJE.cif", "PJE", "C", 5, "PJE")
