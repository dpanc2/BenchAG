#!/usr/bin/env bash
# Usage: bash benchmarks/atac/run_benchmark.sh [output_dir] [smoke_region_count]
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${ALPHAGENOME_API_KEY:?Set ALPHAGENOME_API_KEY before running}"
# Run with Python 3.10 or 3.11 on Linux (TensorFlow 2.15 compatibility).
if [ ! -x .venv/bin/python ]; then bash scripts/setup.sh; fi
.venv/bin/python -m pip install alphagenome==0.7.0 numpy==2.2.6 pandas==2.3.3 scipy==1.15.3 pyBigWig==0.3.29 pyranges==0.1.4 nbclient==0.10.2 nbformat==5.10.4 ipykernel
if [ ! -x data/atac/chrombpnet_env/bin/python ]; then
  mkdir -p data/atac
  .venv/bin/python -m venv data/atac/chrombpnet_env
fi
data/atac/chrombpnet_env/bin/python -m pip install -r benchmarks/atac/requirements.chrombpnet.txt
data/atac/chrombpnet_env/bin/python -m pip install --no-deps chrombpnet==0.1.7
# Reproduce the Colab and generate both model-specific BEDs.
.venv/bin/python -m ipykernel install --prefix "$PWD/.venv" --name benchag
.venv/bin/python - <<'NOTEBOOK'
from pathlib import Path
import nbformat
from nbclient import NotebookClient
root=Path.cwd()
notebook=nbformat.read(root/'benchmarks/atac/prepare_inputs.ipynb',as_version=4)
NotebookClient(notebook,kernel_name='benchag',timeout=600,
               resources={'metadata':{'path':str(root)}}).execute()
nbformat.write(notebook,root/'data/atac/prepare_inputs.executed.ipynb')
NOTEBOOK
# Download the matching IMR-90 signal and model weights.
download() {
  if [ ! -f "$2" ]; then
    curl -L --fail --retry 3 "$1" -o "$2.part"
    mv "$2.part" "$2"
  fi
}
download https://www.encodeproject.org/files/ENCFF650PHA/@@download/ENCFF650PHA.bigWig data/atac/ENCFF650PHA.bigWig
download https://www.encodeproject.org/files/ENCFF113GSV/@@download/ENCFF113GSV.tar.gz data/atac/ENCFF113GSV.tar.gz
printf '%s\n' '6254c06f028aa67f99958c289c8adbed  data/atac/ENCFF650PHA.bigWig' 'f1dd105bcfa444e30b645464f4206fda  data/atac/ENCFF113GSV.tar.gz' | md5sum -c -
mkdir -p data/atac/models data/atac/reference
tar -xzf data/atac/ENCFF113GSV.tar.gz -C data/atac/models --strip-components=2 \
  ./fold_0/model.chrombpnet_nobias.fold_0.ENCSR200OML.h5 \
  ./fold_0/model.bias_scaled.fold_0.ENCSR200OML.h5
if [ ! -f data/atac/models/chrombpnet_recompiled.h5 ]; then
  data/atac/chrombpnet_env/bin/python benchmarks/atac/chrombpnet/reformat_chrombpnet_h5.py \
    -cnb data/atac/models/model.chrombpnet_nobias.fold_0.ENCSR200OML.h5 \
    -bm data/atac/models/model.bias_scaled.fold_0.ENCSR200OML.h5 -o data/atac/models
fi
if [ ! -f data/atac/hg38_chr1_chr3_chr6.fa ]; then
  for chrom in chr1 chr3 chr6; do
    download "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$chrom.fa.gz" "data/atac/reference/$chrom.fa.gz"
  done
  gzip -cd data/atac/reference/chr1.fa.gz data/atac/reference/chr3.fa.gz data/atac/reference/chr6.fa.gz > data/atac/hg38_chr1_chr3_chr6.fa.part
  mv data/atac/hg38_chr1_chr3_chr6.fa.part data/atac/hg38_chr1_chr3_chr6.fa
fi
.venv/bin/python - <<'REFERENCE'
import hashlib
expected={'chr1':'209d95937d1fc25d970390c56f11f79ac9d5ac127ba48c21fbd4e1aac3bb60a5',
          'chr3':'1500ed6a7516eaa5645c56168ff7449f7aadfbd341b94992f63ccdf91b96977c',
          'chr6':'199d76944a1b320cc8e9b3935acd86df00e3bfc3dfa776e6f391f3fc65bf9c6f'}
hashes={}
with open('data/atac/hg38_chr1_chr3_chr6.fa','rb') as f:
    for line in f:
        if line.startswith(b'>'):
            chrom=line[1:].split()[0].decode();hashes[chrom]=hashlib.sha256()
        else:hashes[chrom].update(line.strip().upper().translate(bytes.maketrans(b'RYSWKMBDHV',b'NNNNNNNNNN')))
assert {c:h.hexdigest() for c,h in hashes.items()}==expected,'Unexpected reference sequence'
REFERENCE
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
