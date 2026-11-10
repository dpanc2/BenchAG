#!/usr/bin/env bash
# Usage: bash benchmarks/atac/run_benchmark.sh [output_dir] [smoke_region_count]
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${ALPHAGENOME_API_KEY:?Set ALPHAGENOME_API_KEY before running}"
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
out="${1:-results/atac/run}"
limit=()
if [ -n "${2:-}" ]; then limit=(--limit "$2"); fi
mkdir -p "$out/fold_0" "$out/fold_1" "$out/chrombpnet"
pids=()
for fold in 0 1; do
  .venv/bin/python -u benchmarks/atac/alpha/alpha.py \
    --model-version "FOLD_$fold" --out "$out/fold_$fold" "${limit[@]}" \
    >> "$out/fold_$fold/run.log" 2>&1 &
  pids+=("$!")
done
OMP_NUM_THREADS=4 TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 \
  data/atac/chrombpnet_env/bin/python -u benchmarks/atac/chrombpnet/chrom_bw.py \
  --genome data/atac/hg38_chr1_chr3_chr6.fa --out "$out/chrombpnet" "${limit[@]}" \
  >> "$out/chrombpnet/run.log" 2>&1 &
pids+=("$!")
status=0
for pid in "${pids[@]}"; do if ! wait "$pid"; then status=1; fi; done
if [ "$status" -ne 0 ]; then echo "A run failed; check $out/*/run.log and rerun to resume."; exit 1; fi
.venv/bin/python - "$out" <<'PY'
import json,sys
from pathlib import Path
import pandas as pd
out=Path(sys.argv[1]);rows=[]
for model in ['fold_0','fold_1','chrombpnet']:
    values=json.loads((out/model/'metrics.json').read_text())
    assert values['status']=='complete' and values['n_failed_regions']==0
    rows.append({'run':model,**values})
assert len({r['n_regions'] for r in rows})==1
summary=pd.DataFrame(rows)
summary.to_csv(out/'summary.csv',index=False)
print(summary[['run','n_regions','global_pearson_r_1bp','log_total_count_pearson_r','mean_js_distance']].to_string(index=False))
PY
