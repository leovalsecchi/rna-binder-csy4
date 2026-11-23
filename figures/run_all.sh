#!/usr/bin/env bash
# Rebuild all figure panels. Skips folders named *template*.
set -euo pipefail
cd "$(dirname "$0")"
# The system python3 has none of the deps. Point PY at the venv (see
# environment/README.md) or activate it first.
PY=${PY:-python3}
for d in Fig*/ supp/*/; do
  [ -d "$d" ] || continue
  case "$d" in *template*) continue;; esac
  # every make_*.py, not just make_figure.py: some figures have more than
  # one script (fig03 has make_violins.py, which owns metrics.*).
  shopt -s nullglob
  for s in "$d"make_*.py; do
    echo "== $s"
    "$PY" "$s"
  done
done
