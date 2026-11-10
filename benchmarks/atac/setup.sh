#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
# Run with Python 3.10 or 3.11 on Linux (TensorFlow 2.15 compatibility).
if [ ! -x .venv/bin/python ]; then bash scripts/setup.sh; fi
.venv/bin/python -m pip install alphagenome==0.7.0 numpy==2.2.6 pandas==2.3.3 scipy==1.15.3 pyBigWig==0.3.29
if [ ! -x data/atac/chrombpnet_env/bin/python ]; then
  mkdir -p data/atac
  .venv/bin/python -m venv data/atac/chrombpnet_env
fi
data/atac/chrombpnet_env/bin/python -m pip install -r benchmarks/atac/requirements.chrombpnet.txt
data/atac/chrombpnet_env/bin/python -m pip install --no-deps chrombpnet==0.1.7
.venv/bin/python benchmarks/atac/prepare_inputs.py
