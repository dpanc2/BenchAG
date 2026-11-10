# IMR-90 ATAC-seq benchmark

All code and region tables are included here. Only measurements, model weights,
the reference genome and Python dependencies are downloaded.

## Observed ATAC-seq signal

Ground truth is the IMR-90 ATAC-seq BigWig from
[ENCODE: ENCFF650PHA](https://www.encodeproject.org/files/ENCFF650PHA/).
It is a downloaded input, not included in Git. From the **BenchAG root**:

```bash
mkdir -p data/atac
curl -L --fail --retry 3 \
  https://www.encodeproject.org/files/ENCFF650PHA/@@download/ENCFF650PHA.bigWig \
  -o data/atac/ENCFF650PHA.bigWig
```

The file is 331,977,721 bytes; expected MD5:
`6254c06f028aa67f99958c289c8adbed`.
The benchmark command below also downloads and verifies it automatically if missing.

## Run

Requires Linux, Python 3.10/3.11 and an AlphaGenome API key; no GPU is needed.
From the BenchAG root:

```bash
export ALPHAGENOME_API_KEY='your-key'
bash benchmarks/atac/run_benchmark.sh
```

This prepares the environments/inputs, runs both AG folds and ChromBPNet, and
saves `results/atac/run/summary.csv`. Each model has a log, metrics and cached
profiles. Rerun the same command to resume. For a new run or a two-region test:

```bash
bash benchmarks/atac/run_benchmark.sh results/atac/new_run
bash benchmarks/atac/run_benchmark.sh results/atac/test 2
```

## Code

- `alpha/alpha.py`: AG predictions and metrics; use `--model-version` and `--out`.
- `chrombpnet/chrom_bw.py`: ChromBPNet predictions and evaluation.
- `chrombpnet/chrom_metrics.py`: ChromBPNet metric functions.
- `prepare_inputs.py`: checked downloads, model/reference preparation and fold audit.
- `chrombpnet/reformat_chrombpnet_h5.py`: model-combination helper.
- `data/`: evaluation peaks and genomic folds.
- `requirements.chrombpnet.txt`: CPU environment dependencies.

Based on [the colleague's code, commit 9e6a911](https://github.com/GalinaK413/atac_bench_example/tree/9e6a91157702b80efbc3b43f23c06e21f7bae856).
Prediction/metric calculations are retained. Command-line arguments and cached
profiles replace the old sample-CSV configuration; unused comparison plotting
was removed. The helper comes from ChromBPNet commit `09938fdb` (license beside
it); fold intervals come from Borzoi commit `5c935822`. Execution uses local code.

## Previous full aglab results

All models evaluated 13,812 supplied peak records, using 1 kb midpoint-centered
windows. AG uses 1 Mb input and forward/reverse averaging. ChromBPNet predicts
from summit-centered inputs and averages overlapping outputs; its original
metrics fill uncovered prediction bases with zero.

| Model | Global Pearson | Log-count Pearson | Mean JS distance ↓ |
|---|---:|---:|---:|
| AG FOLD_1, held-out | 0.8664 | 0.8462 | 0.4644 |
| AG FOLD_0, diagnostic | 0.8973 | 0.9077 | 0.4613 |
| ChromBPNet fold0 | 0.7921 | 0.7951 | 0.4722 |

None of the windows overlap FOLD_1 training regions; 13,550 overlap FOLD_0
training regions. Setup records checks in `data/atac/preparation.json`.
Historical draft values differ from these runs.
