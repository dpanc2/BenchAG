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
  --file-mappings data/file_mappings_human.csv --organism human
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

For **mouse**, use `--organism mouse` in all preparation/prediction commands,
mm10 genes, and separate output directories. The ontology/support filenames gain
`_mouse`; AG metadata is `alphagenome_mouse_rna_seq_metadata.csv`.

## Notebooks

- [Prepare the support table from file mappings](notebooks/prepare_ontology_support_from_file_mappings.ipynb).
- [Explain ontology overlap and sample grouping](notebooks/ontology_mapping_explained.ipynb).
- [Check model folds](notebooks/AG_versions.ipynb).

For the existing aglab benchmark, use the [repository quick start](../../README.md).
