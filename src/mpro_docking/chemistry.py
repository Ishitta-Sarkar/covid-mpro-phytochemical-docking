"""Build the co-crystallised 6LU7 inhibitor as one non-covalent molecule.

Current PDB 6LU7 stores inhibitor N3 as chain C (02J-ALA-VAL-LEU-PJE-010),
not as a single residue named N3. A LINK record joins Cys145 SG to PJE C20
(1.79 Å). That bond is omitted here because AutoDock Vina has no covalent
term. Bond orders inside 02J, PJE and 010 come from the RCSB chemical
component definitions in ``data/ccd``. ALA, VAL and LEU use standard amino-acid
connectivity. Consecutive residues are joined by a peptide bond when both
backbone atoms are present; the ester between PJE C and 010 O comes from the
LINK record.

PJE C20-C21 is a double bond in the component dictionary, and the deposited
distance is 1.32 Å, so the ligand keeps that unsaturated warhead. The
crystallographic heavy-atom coordinates are not rebuilt.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from mpro_docking.constants import LIGAND_CHAIN
from mpro_docking.pdbio import AtomRecord, LinkRecord, parse_pdb

BOND_TYPES = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "AROM": Chem.BondType.AROMATIC,
}

# Heavy-atom connectivity for the three standard residues in chain C.
AMINO_ACID_BONDS: dict[str, list[tuple[str, str, str]]] = {
    "ALA": [("N", "CA", "SING"), ("CA", "C", "SING"), ("C", "O", "DOUB"), ("CA", "CB", "SING")],
    "VAL": [
        ("N", "CA", "SING"),
        ("CA", "C", "SING"),
        ("C", "O", "DOUB"),
        ("CA", "CB", "SING"),
        ("CB", "CG1", "SING"),
        ("CB", "CG2", "SING"),
    ],
    "LEU": [
        ("N", "CA", "SING"),
        ("CA", "C", "SING"),
        ("C", "O", "DOUB"),
        ("CA", "CB", "SING"),
        ("CB", "CG", "SING"),
        ("CG", "CD1", "SING"),
        ("CG", "CD2", "SING"),
    ],
}


def parse_ccd_bonds(path: Path) -> list[tuple[str, str, str]]:
    """Return (atom_1, atom_2, order) from a chemical-component CIF."""

    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    bonds: list[tuple[str, str, str]] = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("loop_"):
            i += 1
            continue
        i += 1
        headers: list[str] = []
        while i < len(lines) and lines[i].startswith("_"):
            headers.append(lines[i].split()[0])
            i += 1
        wanted = "_chem_comp_bond.atom_id_1" in headers
        index = {name: n for n, name in enumerate(headers)}
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "loop_", "_")):
            if wanted:
                parts = lines[i].split()
                a = parts[index["_chem_comp_bond.atom_id_1"]]
                b = parts[index["_chem_comp_bond.atom_id_2"]]
                order = parts[index["_chem_comp_bond.value_order"]]
                if not _is_hydrogen_name(a) and not _is_hydrogen_name(b):
                    bonds.append((a, b, order))
            i += 1
    if not bonds:
        raise ValueError(f"no heavy-atom bonds found in {path}")
    return bonds


def _is_hydrogen_name(name: str) -> bool:
    return name.startswith("H")


def _pdb_name(name: str) -> str:
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3s}"


def _load_ccd_dir(ccd_dir: Path) -> dict[str, list[tuple[str, str, str]]]:
    tables: dict[str, list[tuple[str, str, str]]] = {}
    for path in sorted(ccd_dir.glob("*.cif")):
        tables[path.stem] = parse_ccd_bonds(path)
    return tables


def assemble_crystal_ligand(
    pdb_path: Path,
    ccd_dir: Path,
    chain: str = LIGAND_CHAIN,
    relax_h: bool = True,
) -> Chem.Mol:
    """Return the chain-C inhibitor with crystal heavy-atom coordinates and hydrogens."""

    atoms, links = parse_pdb(pdb_path)
    ligand_atoms = [atom for atom in atoms if atom.chain == chain and atom.element != "H"]
    if not ligand_atoms:
        raise ValueError(f"no heavy atoms for chain {chain} in {pdb_path}")
    ccd = _load_ccd_dir(ccd_dir)
    mol = _build_mol(ligand_atoms, links, ccd, chain)
    Chem.SanitizeMol(mol)
    radicals = [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons()]
    if radicals:
        raise ValueError(f"ligand has radical atoms after sanitisation: {radicals}")
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    heavy_xyz = [mol.GetConformer().GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    mol = Chem.AddHs(mol, addCoords=True)
    if relax_h:
        _relax_hydrogens(mol, heavy_xyz)
    mol.SetProp("_Name", "N3")
    return mol


def _build_mol(
    atoms: list[AtomRecord],
    links: list[LinkRecord],
    ccd: dict[str, list[tuple[str, str, str]]],
    chain: str,
) -> Chem.Mol:
    rw = Chem.RWMol()
    index: dict[tuple[str, int, str], int] = {}
    for atom in atoms:
        rd_atom = Chem.Atom(atom.element)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(_pdb_name(atom.name))
        info.SetResidueName(f"{atom.resname:>3s}"[:3])
        info.SetResidueNumber(atom.resseq)
        info.SetChainId(atom.chain)
        rd_atom.SetMonomerInfo(info)
        index[atom.key] = rw.AddAtom(rd_atom)

    seen: set[frozenset[int]] = set()

    def add_bond(key_a: tuple[str, int, str], key_b: tuple[str, int, str], order: str) -> None:
        if key_a not in index or key_b not in index:
            return
        i, j = index[key_a], index[key_b]
        if i == j:
            return
        pair = frozenset((i, j))
        if pair in seen:
            return
        if order not in BOND_TYPES:
            raise ValueError(f"unsupported bond order {order}")
        rw.AddBond(i, j, BOND_TYPES[order])
        seen.add(pair)

    by_residue: dict[tuple[str, int], list[AtomRecord]] = {}
    for atom in atoms:
        by_residue.setdefault((atom.resname, atom.resseq), []).append(atom)

    for (resname, resseq), _residue_atoms in by_residue.items():
        table = AMINO_ACID_BONDS.get(resname, ccd.get(resname))
        if table is None:
            raise ValueError(f"no bond table for residue {resname} {resseq}")
        for name_a, name_b, order in table:
            add_bond((resname, resseq, name_a), (resname, resseq, name_b), order)

    for link in links:
        if link.chain1 != chain or link.chain2 != chain:
            continue
        add_bond(
            (link.resname1, link.resseq1, link.name1),
            (link.resname2, link.resseq2, link.name2),
            "SING",
        )

    ordered = sorted({(atom.resseq, atom.resname) for atom in atoms})
    for (seq_a, name_a), (seq_b, name_b) in zip(ordered, ordered[1:]):
        if seq_b != seq_a + 1:
            continue
        add_bond((name_a, seq_a, "C"), (name_b, seq_b, "N"), "SING")

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom in atoms:
        conf.SetAtomPosition(index[atom.key], atom.xyz)
    mol.AddConformer(conf, assignId=True)
    return mol


def _relax_hydrogens(mol: Chem.Mol, heavy_xyz: list) -> None:
    """Move only hydrogens. Heavy atoms are written back to the crystal positions."""

    try:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
    except Exception:
        ff = None
    if ff is not None:
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() != 1:
                ff.AddFixedPoint(i)
        ff.Minimize(maxIts=200)
    conf = mol.GetConformer()
    heavy_i = 0
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        conf.SetAtomPosition(i, heavy_xyz[heavy_i])
        heavy_i += 1


def residue_atom_key(atom: Chem.Atom) -> tuple[str, int, str]:
    info = atom.GetMonomerInfo()
    if info is None:
        raise ValueError("atom is missing PDB residue information")
    return (info.GetResidueName().strip(), int(info.GetResidueNumber()), info.GetName().strip())
