import math

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from mpro_docking.rmsd import heavy_atom_rmsd, parse_index_map


def _dichloromethane():
    rw = Chem.RWMol()
    carbon = rw.AddAtom(Chem.Atom(6))
    chlorine_a = rw.AddAtom(Chem.Atom(17))
    chlorine_b = rw.AddAtom(Chem.Atom(17))
    rw.AddBond(carbon, chlorine_a, Chem.BondType.SINGLE)
    rw.AddBond(carbon, chlorine_b, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.7, 0.0, 0.0))
    conf.SetAtomPosition(2, Point3D(-1.7, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def test_rmsd_does_not_superimpose():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 1
    assert AllChem.EmbedMolecule(mol, params) == 0
    moved = Chem.Mol(mol)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        point = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(point.x + 3.0, point.y, point.z))
    result = heavy_atom_rmsd(moved, mol)
    assert result.direct == pytest_approx(3.0)
    assert result.symmetry == pytest_approx(3.0)


def test_symmetry_corrects_identical_substituents():
    mol = _dichloromethane()
    swapped = Chem.Mol(mol)
    conf = swapped.GetConformer()
    first = conf.GetAtomPosition(1)
    second = conf.GetAtomPosition(2)
    conf.SetAtomPosition(1, second)
    conf.SetAtomPosition(2, first)
    result = heavy_atom_rmsd(swapped, mol)
    assert result.direct > 1.0
    assert result.symmetry == pytest_approx(0.0)


def test_index_map_parser():
    text = "\n".join(
        [
            "REMARK INDEX MAP 1 5 2 1 3 4",
            "REMARK INDEX MAP 4 2",
            "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 C",
        ]
    )
    assert parse_index_map(text) == {0: 5, 1: 1, 2: 4, 3: 2}


def pytest_approx(value: float):
    import pytest

    return pytest.approx(value, abs=1e-4)


def test_direct_rmsd_matches_manual_shift():
    mol = _dichloromethane()
    moved = Chem.Mol(mol)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        point = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(point.x, point.y + 1.0, point.z))
    result = heavy_atom_rmsd(moved, mol)
    assert math.isclose(result.direct, 1.0, abs_tol=1e-6)
