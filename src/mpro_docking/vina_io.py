"""Write a Vina config, run the binary, and read its scores."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from mpro_docking.constants import (
    CPU,
    ENERGY_RANGE,
    EXHAUSTIVENESS,
    MIN_MODE_RMSD,
    NUM_MODES,
    SEED,
)

_MODE_LINE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$"
)
_VINA_REMARK = re.compile(
    r"REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)"
)
_KCAL_LINE = re.compile(r"(-?\d+\.\d+)\s*\(kcal/mol\)")


@dataclass(frozen=True)
class VinaMode:
    mode: int
    score: float
    rmsd_lb: float
    rmsd_ub: float


def format_vina_config(
    receptor: str,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int = EXHAUSTIVENESS,
    num_modes: int = NUM_MODES,
    energy_range: int = ENERGY_RANGE,
    min_rmsd: int = MIN_MODE_RMSD,
    seed: int = SEED,
) -> str:
    cx, cy, cz = center
    sx, sy, sz = size
    return (
        f"receptor = {receptor}\n"
        f"center_x = {cx:.3f}\n"
        f"center_y = {cy:.3f}\n"
        f"center_z = {cz:.3f}\n"
        f"size_x = {sx:.0f}\n"
        f"size_y = {sy:.0f}\n"
        f"size_z = {sz:.0f}\n"
        f"exhaustiveness = {exhaustiveness}\n"
        f"num_modes = {num_modes}\n"
        f"energy_range = {energy_range}\n"
        f"min_rmsd = {min_rmsd}\n"
        f"seed = {seed}\n"
    )


def parse_vina_modes(text: str) -> list[VinaMode]:
    """Read modes from a Vina log table or from PDBQT ``REMARK VINA RESULT`` lines."""

    remarks = [
        VinaMode(i, float(score), float(lb), float(ub))
        for i, (score, lb, ub) in enumerate(_VINA_REMARK.findall(text), start=1)
    ]
    if remarks:
        return remarks
    modes: list[VinaMode] = []
    for line in text.splitlines():
        match = _MODE_LINE.match(line)
        if match is None:
            continue
        modes.append(
            VinaMode(
                int(match.group(1)),
                float(match.group(2)),
                float(match.group(3)),
                float(match.group(4)),
            )
        )
    return modes


def parse_score_only(text: str) -> float | None:
    for line in text.splitlines():
        if "kcal/mol" not in line:
            continue
        match = _KCAL_LINE.search(line)
        if match:
            return float(match.group(1))
    return None


def torsdof(pdbqt_text: str) -> int:
    match = re.search(r"^TORSDOF\s+(\d+)\s*$", pdbqt_text, flags=re.MULTILINE)
    if match is None:
        raise ValueError("PDBQT file has no TORSDOF record")
    return int(match.group(1))


def run_vina(
    config: Path,
    ligand: Path,
    out: Path | None,
    log: Path,
    cpu: int = CPU,
    extra_args: list[str] | None = None,
    cwd: Path | None = None,
) -> str:
    """Run Vina and write stdout/stderr plus the command to ``log``.

    ``cwd`` must be the repository root when the config uses relative paths.
    """

    log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "vina",
        "--config",
        str(config),
        "--ligand",
        str(ligand),
        "--cpu",
        str(cpu),
    ]
    if extra_args:
        cmd.extend(extra_args)
    elif out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--out", str(out)])
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True, cwd=cwd)
    log.write_text(
        "$ " + " ".join(cmd) + "\n\n" + completed.stdout + "\n" + completed.stderr,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        tail = (completed.stderr or completed.stdout)[-1500:]
        raise RuntimeError(f"vina failed ({completed.returncode}) for {ligand.name}:\n{tail}")
    return completed.stdout + "\n" + completed.stderr
