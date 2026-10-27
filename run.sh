#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TASK="${1:-help}"
if [[ "$#" -gt 0 ]]; then shift; fi
case "$TASK" in
  expression) exec "${PYTHON_BIN:-$ROOT/.venv/bin/python}" "$ROOT/scripts/run_expression.py" "$@" ;;
  mpra|atac|eqtl) echo "$TASK has no contributed implementation yet. See benchmarks/$TASK/README.md" >&2; exit 2 ;;
  *) echo "Usage: ./run.sh expression --config configs/expression.human.yaml [--check]"; exit 2 ;;
esac
