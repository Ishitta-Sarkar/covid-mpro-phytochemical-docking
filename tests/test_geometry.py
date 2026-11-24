from pathlib import Path

from mpro_docking.geometry import centroid, cubic_search_box, distance, is_hydrogen_bond
from mpro_docking.pdbio import parse_pdb

FIXTURE = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "6lu7_chainC.pdb"


def test_box_covers_crystal_ligand():
    atoms, _links = parse_pdb(FIXTURE)
    heavy = [atom.xyz for atom in atoms if atom.chain == "C"]
    center, size = cubic_search_box(heavy)
    assert len(heavy) == 49
    assert center == (-10.729, 12.418, 68.816)
    assert size == (30.0, 30.0, 30.0)
    # Every crystal atom sits inside the rounded box with room to spare.
    for x, y, z in heavy:
        assert abs(x - center[0]) < 15
        assert abs(y - center[1]) < 15
        assert abs(z - center[2]) < 15


def test_centroid_and_distance():
    assert centroid([(0, 0, 0), (2, 0, 0)]) == (1.0, 0.0, 0.0)
    assert distance((0, 0, 0), (3, 4, 0)) == 5


def test_hydrogen_bond_geometry():
    donor = (0.0, 0.0, 0.0)
    hydrogen = (1.0, 0.0, 0.0)
    acceptor = (3.0, 0.0, 0.0)
    ok, dist, ang = is_hydrogen_bond(donor, hydrogen, acceptor)
    assert ok
    assert round(dist, 3) == 2.0
    assert round(ang, 1) == 180.0

    too_far = (4.0, 0.0, 0.0)
    ok, dist, _ang = is_hydrogen_bond(donor, hydrogen, too_far)
    assert not ok
    assert dist > 2.5

    bent = (1.0, 2.0, 0.0)
    ok, _dist, ang = is_hydrogen_bond(donor, hydrogen, bent)
    assert not ok
    assert ang < 120
