# IMR-90 ATAC-seq benchmark

All benchmark code and region tables are included here. Setup downloads only
measurements, model weights, the reference genome and Python dependencies.
It does not clone repositories or download Python scripts from GitHub.

## Files

- `alpha/alpha.py`: colleague's AlphaGenome prediction and metric functions.
- `chrombpnet/chrom_bw.py`: colleague's ChromBPNet prediction code.
- `chrombpnet/chrom_metrics.py`: colleague's metric code.
- `chrombpnet/reformat_chrombpnet_h5.py`: official model-combination helper.
- `data/`: evaluation peaks and genomic fold intervals.
- `run_alphagenome.py`, `run_chrombpnet.py`: explicit paths/folds and resumable caches.
- `prepare_inputs.py`, `setup.sh`: download/check inputs and prepare environments.

The three colleague scripts and evaluation BED are copied unchanged from
[atac_bench_example, commit 9e6a911](https://github.com/GalinaK413/atac_bench_example/tree/9e6a91157702b80efbc3b43f23c06e21f7bae856).
The helper is copied from ChromBPNet commit `09938fdb`; its license is in
`chrombpnet/LICENSE.reformat_helper`. Fold intervals come from Borzoi commit
`5c935822`. These links record provenance; execution uses the local files.
The wrapper omits the original unused torch import when loading AlphaGenome code.

## Run

Requires Linux, Python 3.10 or 3.11, Git for initial BenchAG checkout, internet
access and an AlphaGenome API key. No GPU is needed. From the BenchAG root:

```bash
source .venv/bin/activate
export ALPHAGENOME_API_KEY='your-key'
bash benchmarks/atac/run_benchmark.sh
```

If `.venv` does not exist, omit the activation line; the command creates it.
Setup installs dependencies, downloads/checks public inputs and prepares the
ChromBPNet model. Results go to **`results/atac/run/summary.csv`**, with a log,
metrics and cached predictions in each model folder. Previous runs are preserved.
Rerun with the same output folder to resume; choose a new folder for fresh predictions:

```bash
bash benchmarks/atac/run_benchmark.sh results/atac/new_run
```

For two regions only, add `2` as the final argument. To prepare inputs without
running predictions: `bash benchmarks/atac/setup.sh`.

## Evaluation and previous results

All models use the same 13,812 supplied peak records. AG uses a 1 Mb context and
forward/reverse averaging. Metrics use midpoint-centered 1 kb windows.
ChromBPNet uses summit-centered inputs and averages overlapping outputs. Its
original metric code fills uncovered prediction bases with zero.

| Model | Global Pearson | Log-count Pearson | Mean JS distance ↓ |
|---|---:|---:|---:|
| AG FOLD_1, held-out | 0.8664 | 0.8462 | 0.4644 |
| AG FOLD_0, diagnostic | 0.8973 | 0.9077 | 0.4613 |
| ChromBPNet fold0 | 0.7921 | 0.7951 | 0.4722 |

These full runs were completed on aglab. Use FOLD_1 for held-out evaluation:
none of its evaluation windows overlap training regions; 13,550 overlap FOLD_0
training regions. Input checks and model-input reference checksums are saved in
`data/atac/preparation.json`. The public reference is checked after the source
encoder's ambiguous-base handling. The historical draft values differ from
these runs; fresh API results must be compared rather than assumed identical.
