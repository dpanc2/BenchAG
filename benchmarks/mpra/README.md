# CAGI5 MPRA benchmark

Reproduce [tsoiegor/alphabench](https://github.com/tsoiegor/alphabench) at commit
`58ed6e7a8d52cf88fdfea4125f5db21f723c5623`, using the Enformer and Borzoi
ontology mappings. Both run AlphaGenome FOLD_1: these are mapping profiles,
not predictions from Enformer or Borzoi.

The supplied config keeps the original procedure: Tags > 10, source Position + 1,
1 Mb input, DNASE CenterMaskScorer (width 501, DIFF_SUM), and mean over matching
tracks. It evaluates signed Pearson and Spearman correlations against reporter
activity, overall and per (Element, Cell_Type) assay. Thus this is a DNase-based
proxy for MPRA variant effects, not a direct reporter-expression output.

Use a Linux GPU environment with Python 3.11+, CUDA-enabled JAX and AlphaGenome
Research, and access to the FOLD_1 weights. The expression API environment alone
is insufficient. Install the benchmark:

```bash
python -m pip install 'git+https://github.com/tsoiegor/alphabench.git@58ed6e7a8d52cf88fdfea4125f5db21f723c5623'
alphabench download --config benchmarks/mpra/cagi5.yaml
alphabench process --config benchmarks/mpra/cagi5.yaml
```

After resolving the reference mismatch below, run from the BenchAG root:

```bash
CUDA_VISIBLE_DEVICES=0 nohup alphabench run --config benchmarks/mpra/cagi5.yaml \
  > results/mpra.log 2>&1 < /dev/null &
```

Outputs go to `results/mpra/<run>/`, including predictions.tsv, predictions.parquet,
raw scores, metrics.json and manifests for each mapping. Choose an available GPU;
with two visible GPUs, the profiles can run in parallel. Resume with
`--resume results/mpra/<run>`.

## Reproduction status

Local data preparation reproduced both profiles: 44,647 input rows, 29,035 selected
rows, 20,421 unique variants, 14 loci and 22 assay groups. The input MD5 matches
`eacbb118660a3bd39575d88b1a03f7a7`.

Independent UCSC hg38 checks found one REF mismatch: source row_id 27570,
SORT1-flip, chr1:109274650 (zero-based), T deletion; hg38 has G. The upstream
validator currently prints mismatches and continues. No allele correction or
exclusion has been applied; resolve and record this before the paper run.

Preparation tables and the reference audit are saved locally under
`results/mpra/preflight/` and are not prediction results. GPU inference has not
started because aglab SSH is unavailable. No MPRA scores have been added to the paper.
