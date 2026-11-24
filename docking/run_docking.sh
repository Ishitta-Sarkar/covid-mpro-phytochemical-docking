#!/usr/bin/env bash
# Fetch PDB 6LU7 and the PubChem ligands, prepare them, and run AutoDock Vina.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"

if ! command -v vina >/dev/null 2>&1; then
  echo "vina was not found on PATH." >&2
  echo "Create the environment with: conda env create -f environment.yml" >&2
  echo "Then activate it: conda activate mpro-docking" >&2
  exit 1
fi
if ! command -v mk_prepare_receptor.py >/dev/null 2>&1; then
  echo "mk_prepare_receptor.py was not found on PATH. Activate the mpro-docking environment." >&2
  exit 1
fi

exec python -m mpro_docking "$@"
