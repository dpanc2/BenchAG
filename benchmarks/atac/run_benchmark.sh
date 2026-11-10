#!/usr/bin/env bash
# Usage: bash benchmarks/atac/run_benchmark.sh [output_dir] [smoke_region_count]
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${ALPHAGENOME_API_KEY:?Set ALPHAGENOME_API_KEY before running}"
if [ ! -f data/atac/fold1_test_regions.bed ]; then
  echo "Run benchmarks/atac/prepare_inputs.ipynb first."
  exit 1
fi
# Run with Python 3.10 or 3.11 on Linux (TensorFlow 2.15 compatibility).
if [ ! -x .venv/bin/python ]; then bash scripts/setup.sh; fi
.venv/bin/python -m pip install alphagenome==0.7.0 numpy==2.2.6 pandas==2.3.3 scipy==1.15.3 pyBigWig==0.3.29
if [ ! -x data/atac/chrombpnet_env/bin/python ]; then
  mkdir -p data/atac
  .venv/bin/python -m venv data/atac/chrombpnet_env
fi
data/atac/chrombpnet_env/bin/python -m pip install -r benchmarks/atac/requirements.chrombpnet.txt
data/atac/chrombpnet_env/bin/python -m pip install --no-deps chrombpnet==0.1.7
# Download the matching IMR-90 signal and model weights.
download() {
  if [ ! -f "$2" ]; then
    curl -L --fail --retry 3 "$1" -o "$2.part"
    mv "$2.part" "$2"
  fi
}
download https://www.encodeproject.org/files/ENCFF650PHA/@@download/ENCFF650PHA.bigWig data/atac/ENCFF650PHA.bigWig
download https://www.encodeproject.org/files/ENCFF113GSV/@@download/ENCFF113GSV.tar.gz data/atac/ENCFF113GSV.tar.gz
mkdir -p data/atac/models data/genomes
tar -xzf data/atac/ENCFF113GSV.tar.gz -C data/atac/models --strip-components=2 \
  ./fold_0/model.chrombpnet_nobias.fold_0.ENCSR200OML.h5 \
  ./fold_0/model.bias_scaled.fold_0.ENCSR200OML.h5
if [ ! -f data/atac/models/chrombpnet_recompiled.h5 ]; then
  data/atac/chrombpnet_env/bin/python benchmarks/atac/chrombpnet/reformat_chrombpnet_h5.py \
    -cnb data/atac/models/model.chrombpnet_nobias.fold_0.ENCSR200OML.h5 \
    -bm data/atac/models/model.bias_scaled.fold_0.ENCSR200OML.h5 -o data/atac/models
fi
# One full hg38 reference, shared by the human benchmarks.
if [ ! -f data/genomes/hg38.fa ]; then
  download https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz data/genomes/hg38.fa.gz
  gzip -dc data/genomes/hg38.fa.gz > data/genomes/hg38.fa.part
  mv data/genomes/hg38.fa.part data/genomes/hg38.fa
fi
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
  --genome data/genomes/hg38.fa --out "$out/chrombpnet" "${limit[@]}" \
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
    rows.append({'run':model,**values})
summary=pd.DataFrame(rows)
summary.to_csv(out/'summary.csv',index=False)
print(summary[['run','n_regions','global_pearson_r_1bp','log_total_count_pearson_r','mean_js_distance']].to_string(index=False))
PY
