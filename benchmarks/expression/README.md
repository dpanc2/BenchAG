# Expression benchmark

Run commands from the **BenchAG repository root** with the environment activated.
First-time setup: `bash scripts/setup.sh`. Set `ALPHAGENOME_API_KEY` before making
AlphaGenome requests.

## 1. Match your samples to AlphaGenome

Start with your original `file_mappings.csv`: it needs an `id` column containing
ENCODE accessions. Keep `original_id` if available; it provides a fallback accession.
For mixed-species files, `genome` or `assembly` selects the requested organism.

```bash
python benchmarks/expression/scripts/build_alphagenome_ontology_map.py \
  --file-mappings data/file_mappings_borzoi_and_v1_human.csv --organism human
```

This looks up the ontology and assay in ENCODE, then compares them with AG.
It writes these files under `benchmarks/expression/data/`:

- `qnorm_alphagenome_ontology_map.csv`: ENCODE metadata for each sample.
- `qnorm_alphagenome_support_summary.csv`: whether AG has that ontology and assay.

Check **`alphagenome_ontology_assay_supported`** for an exact match.
`alphagenome_supported` checks the ontology alone. Unresolved samples stay in the
table with their errors. ENCODE responses are cached; rerunning reuses them.
AG metadata is read from `alphagenome_human_rna_seq_metadata.csv`, or exported if
absent. `--metadata` selects a different snapshot. The historical `qnorm` filenames
do not mean this step normalizes expression.

Create the JSON mappings used below:

```bash
python scripts/prepare_expression_mapping.py \
  --ontology-map benchmarks/expression/data/qnorm_alphagenome_ontology_map.csv \
  --metadata benchmarks/expression/data/alphagenome_human_rna_seq_metadata.csv \
  --out inputs/human_mapping
```

## 2. Predict expression

Supply `inputs/genes.tsv` with columns `gene_id`, `chrom`, `TSS`, `gene_strand`.
Use hg38 coordinates, a zero-based TSS, and strand `+` or `-`.

```bash
python scripts/predict_expression.py \
  --genes inputs/genes.tsv \
  --mapping inputs/human_mapping/track_to_first_id.json \
  --organism human --window 128kb --out results/my_run
```

Windows: `16kb`, `128kb`, or `1mb`. The default model is `FOLD_1`.
Results are in `results/my_run/predictions.csv`, with genes as rows and
ontology–assay names as columns. Missing-annotation genes are recorded and skipped.
Rerun the same command to resume. An all-NaN column means no usable output under
the gene-strand selection; it does not mean zero expression.

## 3. Compare with measured expression

Supply `inputs/truth.csv` with gene IDs as rows and experimental sample IDs as
columns. The transposed layout is also accepted.

```bash
python scripts/compare_expression.py \
  --predictions results/my_run/predictions.csv --truth inputs/truth.csv \
  --mapping inputs/human_mapping/track_to_all_ids.json \
  --out results/my_run/comparison
```

The command applies `log2(x + 1)` to measurements, then averages samples sharing
an ontology and assay. Use `--truth-transform none` for already transformed truth.
It reports gene and cell Pearson correlations in `comparison/summary.json`.

The included `data/file_mappings_borzoi_all.csv` contains both human and mouse samples.
For **mouse**, use that file and `--organism mouse` in all preparation/prediction commands,
mm10 genes, and separate output directories. The ontology/support filenames gain
`_mouse`; AG metadata is `alphagenome_mouse_rna_seq_metadata.csv`.

## Prepare ground-truth TPM values

The preparation code is adapted from Benjamin's
[ExpressionPredict scripts](https://github.com/minjaf/ExpressionPredict/tree/6c10d3b7f2022488f40c6d3c52ddca338756cea5/src).
It downloads the ENCODE `encode_tsv` accessions recorded in file mappings, builds
a gene-by-sample TPM matrix, and quantile-normalizes each species separately.
BAM files are not needed.

Use a separate environment for the original normalization dependency (NumPy < 2):

```bash
python -m venv --system-site-packages data/gt_env
source data/gt_env/bin/activate
pip install 'qnorm @ git+https://github.com/Maarten-vd-Sande/qnorm.git@9789ed3'
```

```bash
python benchmarks/expression/scripts/prepare_ground_truth.py \
  --file-mappings data/file_mappings_borzoi_and_v1_human.csv \
  --organism human \
  --annotations inputs/GRCh38_genes.csv \
  --intervals-prefix inputs/human. \
  --out results/ground_truth/human
```

Supply the original annotation CSV (`gene_id`, `chromosome`, `TSS`, `TES`) and
six tab-separated interval tables: `{train,valid,test}.{forward,reverse}.csv`,
with `gene_id`, `chrom`, `TSS`, `TES`. Gene IDs and coordinates must match.
For mouse, use `--organism mouse`, mm10 annotations/intervals, and its own output
folder; `data/file_mappings_borzoi_all.csv` contains mouse samples.

Outputs: `human_tpms.csv`, `human_qnorm_tpms.csv`, sample list and a manifest
with ENCODE accessions and input checksums. Downloads are cached on reruns.
Use either matrix as `--truth` in the comparison command above.

**Historical qnorm values require the original full normalization pool.**
The supplied mapping files are useful starting points, but do not include all
additional samples listed in Benjamin's original commands. Pass all original
mapping files together to reproduce that pool. Normalize before AG matching;
changing the sample or gene pool changes qnorm values. The script preserves the
original annotation/interval intersection and training-gene filtering rule;
all validation/test genes surviving the intersection are retained. Like the original downloader, it removes the first quantification row.

## Notebooks

- [Prepare the support table from file mappings](notebooks/prepare_ontology_support_from_file_mappings.ipynb).
- [Explain ontology overlap and sample grouping](notebooks/ontology_mapping_explained.ipynb).
- [Model folds and gene TSS verification](notebooks/AG_versions.ipynb).

For the existing aglab benchmark, use the [repository quick start](../../README.md).
