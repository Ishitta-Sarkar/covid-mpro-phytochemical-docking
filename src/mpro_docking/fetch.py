"""Download the receptor and ligands from RCSB PDB and PubChem."""

from __future__ import annotations

import hashlib
import urllib.error
import urllib.request
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from mpro_docking.constants import PDB_URL, PUBCHEM_2D_URL, PUBCHEM_3D_URL, SEED, USER_AGENT


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def download(url: str, dest: Path, refresh: bool = False) -> Path:
    if dest.exists() and not refresh and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise urllib.error.HTTPError(
            exc.url, exc.code, f"{exc.reason} ({url})", exc.hdrs, exc.fp
        ) from exc
    if not payload:
        raise RuntimeError(f"empty response from {url}")
    dest.write_bytes(payload)
    return dest


def fetch_pdb(pdb_id: str, dest: Path, refresh: bool = False) -> Path:
    return download(PDB_URL.format(pdb_id=pdb_id), dest, refresh=refresh)


def fetch_pubchem_sdf(cid: str, dest: Path, refresh: bool = False) -> tuple[Path, str]:
    """Download a 3D SDF, falling back to a 2D SDF if PubChem has no 3D record."""

    try:
        download(PUBCHEM_3D_URL.format(cid=cid), dest, refresh=refresh)
        source = "pubchem-3d"
    except urllib.error.HTTPError:
        download(PUBCHEM_2D_URL.format(cid=cid), dest, refresh=True)
        source = "pubchem-2d"
    return dest, source


def load_ligand_sdf(path: Path, source: str, seed: int = SEED) -> tuple[Chem.Mol, str]:
    """Load one ligand, embedding a conformer only when the file has no coordinates."""

    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    mol = next((item for item in supplier if item is not None), None)
    if mol is None:
        raise ValueError(f"RDKit could not read {path}")
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise RuntimeError(f"ETKDGv3 failed to embed {path}")
        source = f"{source}+rdkit-etkdg-v3-seed-{seed}"
    elif not any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        mol = Chem.AddHs(mol, addCoords=True)
    Chem.SanitizeMol(mol)
    return mol, source
