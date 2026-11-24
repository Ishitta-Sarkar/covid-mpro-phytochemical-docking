"""Repository paths. Entry points resolve these from this file, not from cwd."""

from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def data_dir() -> Path:
    return repo_root() / "data"


def results_dir() -> Path:
    return repo_root() / "results"
