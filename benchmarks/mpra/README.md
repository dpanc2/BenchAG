# CAGI5 MPRA benchmark

Uses [tsoiegor/alphabench](https://github.com/tsoiegor/alphabench), pinned to
`58ed6e7a8d52cf88fdfea4125f5db21f723c5623`. The Enformer and Borzoi profiles
are ontology mappings; both run AlphaGenome ALL_FOLDS.

Settings: Tags > 10, source Position + 1, 1 Mb input, DNASE CenterMaskScorer
(width 501, DIFF_SUM), mean over matching tracks. Metrics are signed Pearson
and Spearman correlations with reporter activity, overall and per
(Element, Cell_Type). DNase variant effects serve as a proxy for MPRA effects.

## Setup

Use Python 3.11+, CUDA-enabled JAX, AlphaGenome Research, and access to ALL_FOLDS
weights. Clone the pinned package into an environment with those dependencies:

```bash
git clone https://github.com/tsoiegor/alphabench.git data/alphabench
cd data/alphabench
git checkout 58ed6e7a8d52cf88fdfea4125f5db21f723c5623
git apply ../../benchmarks/mpra/alphabench-splice-reference.patch
python -m pip install -e .
cd ../..
```

The patch supplies standard splice-site annotations needed by AlphaGenome
Research 0.3.0. Without them, variant inference fails at `jnp.asarray(None)`.
It does not change the DNase scorer or aggregation.

## Shared genome

The config uses `data/genomes/hg38.fa`, the same full hg38 reference as ATAC.
Place the reference there before running MPRA. On aglab it is already provided;
on a new machine you can download it from UCSC:

```bash
mkdir -p data/genomes
curl -L --fail --retry 3 \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  -o data/genomes/hg38.fa.gz
gzip -dc data/genomes/hg38.fa.gz > data/genomes/hg38.fa
```

## Run

```bash
alphabench download --config benchmarks/mpra/cagi5.yaml
alphabench process --config benchmarks/mpra/cagi5.yaml
CUDA_VISIBLE_DEVICES=0 nohup alphabench run --config benchmarks/mpra/cagi5.yaml \
  > results/mpra.log 2>&1 < /dev/null &
```

Resolve the reference mismatch below before a paper run. Output folders contain
`predictions.tsv`, `predictions.parquet`, raw scores, `metrics.json`, manifests,
and saved batches. Resume with `--resume results/mpra/<run>` and the same config.

## Reproduction on anogena

Workspace: `/home/jovyan/dpanc/BenchAG`. Environment:
`/home/jovyan/miniconda3/envs/alphagenome/bin/python`.
The previous FOLD_1 run is complete in
`results/mpra/20261103T214431Z-04c334bd/`.

The ALL_FOLDS rerun is queued until GPU 0 is free. It uses the same filtered
input, two scoring threads, batches of 16 and sequential mapping profiles.

- Config: `results/mpra/all_folds/cagi5.yaml`.
- Queue/status: `results/mpra/all_folds/current_run.json`.
- Log after inference starts: `results/mpra/all_folds/full.log`.
- Outputs: `results/mpra/all_folds/runs/`.
- Comparison after completion: `results/mpra/all_folds/comparison_with_fold1.csv`.

A small ALL_FOLDS test runs automatically before the full dataset. Existing
training and GENA inference jobs are not interrupted.

The original input contains 44,647 rows; its MD5 is
`eacbb118660a3bd39575d88b1a03f7a7`. Strict checks against hg38 found one bad REF:
source row_id 27570, SORT1-flip, chr1:109274650 zero-based, T deletion; hg38 has G.
The upstream validator only prints this mismatch. We exclude that row, retaining
its original ID as `source_row_id`, and record the exclusion in
`results/mpra/reference_filter.json`. No allele is changed.

After exclusion and the original filters, each profile has **29,034 rows,
20,420 unique variants, 14 loci and 22 assay groups**. To recreate the filtered
input from the downloaded original table:

```python
import pandas as pd
f = pd.read_csv('data/mpra/Kircher_GRCh38_ALL.tsv', sep='\t')
f.insert(0, 'source_row_id', range(len(f)))
f = f[f.source_row_id != 27570]
f.to_csv('data/mpra/Kircher_GRCh38_REF_matched.tsv', sep='\t', index=False)
```

In a copy of the config, set `variants_url` to the filtered file's absolute
`file://` URL and `variants_checksum` to its `sha256:<digest>`. The saved run config
records both. The anogena run uses its validated existing hg38 FASTA, two scoring
threads, batches of 16, and sequential profiles. Scoring settings are unchanged.
The checkpoint revision, package versions, source/patch hashes and reference
checksum are recorded in `results/mpra/reproduction_environment.json`.

## Comparison with the paper

ALL_FOLDS matches the paper's model choice, but the current pipeline is not an
exact reproduction of every [paper setting](https://www.biorxiv.org/content/10.1101/2025.06.25.661532v1.full).
It uses hg38 and a 501 bp DNase ALT-minus-REF score for both mapping profiles.
The linked preprint methods describe hg19, a 512 bp difference for the Enformer
strategy and a 4,000 bp log-fold-change for the Borzoi strategy, plus specific
locus exclusions. Mapping names in this config change tissues only.
The [final paper](https://www.nature.com/articles/s41586-025-10014-0) reports
zero-shot cell-type-matched DNase Pearson r ≈ 0.57
(Fig. 5j, averaged across loci). This is a comparison target, not a result of our
queued run. Report pooled, mean-per-assay and mean-per-locus correlations separately;
changing the checkpoint alone does not resolve those protocol differences.
