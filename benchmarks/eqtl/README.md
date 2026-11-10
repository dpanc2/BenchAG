# eQTL datasets

Source: [calengil/AlphaGenome-eQTL](https://github.com/calengil/AlphaGenome-eQTL),
commit `79c4cdf144a2ea21ec249b6e0d61902ce73e6ebb`.
The original datasets and scripts are used unchanged.

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
