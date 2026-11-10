# eQTL datasets

Source: [calengil/AlphaGenome-eQTL](https://github.com/calengil/AlphaGenome-eQTL),
commit `79c4cdf144a2ea21ec249b6e0d61902ce73e6ebb`.
The supplied datasets and scripts are retained in `data/eqtl/source`; the wrappers below rebuild datasets and add resumable prediction runs.

From the BenchAG root:

```bash
mkdir -p data/eqtl
git clone https://github.com/calengil/AlphaGenome-eQTL.git data/eqtl/source
git -C data/eqtl/source checkout --detach 79c4cdf144a2ea21ec249b6e0d61902ce73e6ebb
```

| File in `data/eqtl/source/sets/` | Variant–gene–tissue pairs | Unique variants |
|---|---:|---:|
| `snv_eqtl.tsv` | 17,675 | 6,208 |
| `indel_credible_sets_split.tsv` | 2,673 | 829 |
| `indel_credible_sets_merged.tsv` | 2,645 | 829 |

All three cover 49 GTEx tissues. `tissue_mapping.tsv` maps dataset IDs to AG
tissue labels. The indel versions treat the cerebellum datasets separately or
merge them. These are supplied prepared datasets, not a reconstruction from
raw eQTL Catalogue files.

## Check saved predictions

```bash
python -m venv --system-site-packages data/eqtl/env
data/eqtl/env/bin/pip install pandas pyarrow scipy scikit-learn
data/eqtl/env/bin/python data/eqtl/source/code/eqtl_metrics.py \
  data/eqtl/source/feather/snv_eqtl_ALL_FOLDS_sign.feather \
  data/eqtl/source/feather/snv_eqtl_ALL_FOLDS_coeff.feather
```

Replace the filename prefix with another dataset/model to evaluate it.
This uses the colleague's saved predictions; it does not make new AG requests.
Metrics are calculated per tissue, then weighted by the number of valid rows.

On aglab, the dataset audit is in `results/eqtl/dataset_audit.json` and metrics
for all supplied ALL_FOLDS/FOLD_1 predictions are in
`results/eqtl/saved_prediction_metrics.csv`.

## Rebuild datasets and generate fresh predictions

With the BenchAG environment activated (requires `pyarrow`):

```bash
python benchmarks/eqtl/prepare_datasets.py
python benchmarks/eqtl/run_predictions.py \
  --variants data/eqtl/rebuilt/indel_credible_sets_split.tsv \
  --mapping data/eqtl/rebuilt/tissue_mapping.tsv \
  --model ALL_FOLDS --out results/eqtl/fresh/ALL_FOLDS/indel_credible_sets_split
```

Set `ALPHAGENOME_API_KEY` for prediction. Other datasets: `snv_eqtl` and
`indel_credible_sets_merged`; other model: `FOLD_1`. Rerun the same command to
resume from cached API scores. Use a separate output folder for each dataset/model.

Preparation reconstructs SNV targets from the
[official evaluation data](https://storage.googleapis.com/alphagenome/evals/eqtl_variant_borzoi_coefficient_human_predictions.feather).
Indels are downloaded from raw GTEx gene-expression credible sets in the
[eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/), then filtered to non-SNV
variants with PIP ≥ 0.9 and genes in the official GENCODE v46 annotation.
Duplicate variant/gene rows within a dataset are removed. The merged version
keeps the smallest-p-value record per variant/gene/tissue label. Both variants
of cerebellum remain separate in the split version.

`data/eqtl/rebuilt/preparation_manifest.json` records raw-file URLs/checksums and
checks pair identities and effect sizes against the colleague's supplied tables.
Predictions run the original `make_preds.py` scoring, aggregation and export code;
the API/tidy block is wrapped to cache the required gene/tissue scores. Missing
scores retain the original zero fallback. Fresh cache files never load the
colleague's saved predictions. These files are used only for later comparison.
