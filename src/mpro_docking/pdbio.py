"""Minimal PDB reader for ATOM, HETATM and LINK records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AtomRecord:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    x: float
    y: float
    z: float
    element: str
    line: str

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.resname, self.resseq, self.name)


@dataclass(frozen=True)
class LinkRecord:
    name1: str
    resname1: str
    chain1: str
    resseq1: int
    name2: str
    resname2: str
    chain2: str
    resseq2: int


def _field(line: str, start: int, end: int) -> str:
    """Slice a 1-based inclusive PDB column range."""

    return line[start - 1 : end].strip()


def _element(line: str, name: str) -> str:
    if len(line) >= 78 and line[76:78].strip():
        return line[76:78].strip().title()
    letters = "".join(ch for ch in name if ch.isalpha())
    if not letters:
        return ""
    if len(letters) > 1 and letters[1].islower():
        return letters[:2].title()
    return letters[0].upper()


def parse_atom(line: str) -> AtomRecord:
    name = _field(line, 13, 16)
    return AtomRecord(
        record=line[0:6].strip(),
        serial=int(_field(line, 7, 11) or "0"),
        name=name,
        altloc=_field(line, 17, 17),
        resname=_field(line, 18, 20),
        chain=_field(line, 22, 22),
        resseq=int(_field(line, 23, 26) or "0"),
        x=float(_field(line, 31, 38)),
        y=float(_field(line, 39, 46)),
        z=float(_field(line, 47, 54)),
        element=_element(line, name),
        line=line.rstrip("\n"),
    )


def parse_link(line: str) -> LinkRecord:
    return LinkRecord(
        name1=_field(line, 13, 16),
        resname1=_field(line, 18, 20),
        chain1=_field(line, 22, 22),
        resseq1=int(_field(line, 23, 26) or "0"),
        name2=_field(line, 43, 46),
        resname2=_field(line, 48, 50),
        chain2=_field(line, 52, 52),
        resseq2=int(_field(line, 53, 56) or "0"),
    )


def parse_pdb(path: Path) -> tuple[list[AtomRecord], list[LinkRecord]]:
    """Parse the first model. Alternate locations other than A are dropped."""

    atoms: list[AtomRecord] = []
    links: list[LinkRecord] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_atom(line)
            if atom.altloc not in ("", "A"):
                continue
            atoms.append(atom)
        elif line.startswith("LINK"):
            links.append(parse_link(line))
    return atoms, links


def write_atom_lines(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(line.rstrip() for line in lines + ["TER", "END"]) + "\n"
    path.write_text(text, encoding="utf-8")
