# IMR-90 ATAC-seq benchmark

All three full runs were completed on aglab: AlphaGenome FOLD_0, FOLD_1 and
ChromBPNet fold0, each evaluated on **13,812 peak records**.

## Run from a fresh checkout

Requires Linux, Python 3.10 or 3.11, Git, internet access and an AlphaGenome API
key. No GPU is needed. Run from the BenchAG root:

```bash
export ALPHAGENOME_API_KEY='your-key'
bash benchmarks/atac/run_benchmark.sh
```

The command installs the environments, downloads and checks inputs, prepares
the model/reference, runs all three models, and writes `results/atac/summary.csv`.
Each model has its own `run.log`, `metrics.json` and cached predictions.
Rerun the same command to resume. For a quick test on two regions:

```bash
bash benchmarks/atac/run_benchmark.sh results/atac/smoke 2
```

To prepare inputs without inference: `bash benchmarks/atac/setup.sh`.
The first run downloads about 1.3 GB of inputs plus Python dependencies; allow
several GB of disk space for references, environments and predictions.

## What it reproduces

Prediction and metric functions come from
[the colleague's repository](https://github.com/GalinaK413/atac_bench_example/tree/9e6a91157702b80efbc3b43f23c06e21f7bae856).
Wrappers add explicit folds, input manifests and resumable prediction caches.
The unused upstream torch import is omitted.

Inputs are ENCODE `ENCFF650PHA` (observed signal), `ENCFF113GSV` (models), and
`ENCFF802KJC` (original regions), with checked MD5 values. The pinned source's
`alpha/fold1_test_regions.bed` supplies the evaluation peaks; setup verifies that
they belong to ChromBPNet's original fold0 test peaks. Duplicate evaluation
windows are retained. The public UCSC hg38 chr1/chr3/chr6 sequences are verified
against our completed-run reference after mapping ambiguous letters to `N`,
as the source one-hot encoder does. The original CPU model is rebuilt
from the no-bias and scaled-bias weights using the pinned official helper.
Input checks and fold-overlap counts are saved in `data/atac/preparation.json`.

AlphaGenome uses a 1 Mb context, IMR-90 (`EFO:0001196`), and forward/reverse
averaging. Evaluation uses 1 kb windows at peak midpoints. ChromBPNet predicts
from summit-centered 2,114 bp inputs and averages overlapping outputs into a
BigWig. Its original metric code fills uncovered prediction bases with zero;
our full run has 239,797 such base slots. These original choices are retained.

## Completed aglab results

| Model | Global Pearson | Log-count Pearson | Mean JS distance ↓ |
|---|---:|---:|---:|
| AG FOLD_1, held-out | 0.8664 | 0.8462 | 0.4644 |
| AG FOLD_0, diagnostic | 0.8973 | 0.9077 | 0.4613 |
| ChromBPNet fold0 | 0.7921 | 0.7951 | 0.4722 |

Use **FOLD_1** for the held-out comparison: none of the evaluation windows overlap
its training regions. For FOLD_0, 13,550 overlap training regions, so its result
is diagnostic. These are our reproduced-run results; they differ from the
historical draft values. Fresh API predictions may vary; compare measured
results rather than assuming bit-for-bit equality.

Individual commands remain available in `run_alphagenome.py --help` and
`run_chrombpnet.py --help`. ChromBPNet also saves its predicted BigWig.
