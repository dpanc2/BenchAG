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

## Prepare your own ground truth

Two scripts, based on Benjamin's [ExpressionPredict code](https://github.com/minjaf/ExpressionPredict/tree/6c10d3b7f2022488f40c6d3c52ddca338756cea5/src):

```bash
python benchmarks/expression/scripts/download_tpms.py \
  --file-mappings data/file_mappings_borzoi_and_v1_human.csv \
  --organism human --out results/ground_truth/human
```

This downloads the `encode_tsv` accessions from file mappings, caches the files,
and writes `human_tpms.csv` (genes × samples). It keeps Ensembl genes present in
every sample. No BAMs, gene annotation or fold tables are needed.

For normalization, use a separate environment for the original qnorm dependency:

```bash
python -m venv --system-site-packages data/gt_env
source data/gt_env/bin/activate
pip install 'qnorm @ git+https://github.com/Maarten-vd-Sande/qnorm.git@9789ed3'
python benchmarks/expression/scripts/qnorm_tpms.py \
  --tpms results/ground_truth/human/human_tpms.csv \
  --out results/ground_truth/human/human_qnorm_tpms.csv
```

You can also normalize your own CSV: `gene_id` followed by numeric sample columns,
with no missing values. Normalization uses raw TPMs, before the comparison step's
log transform. It handles ties and keeps all input genes.

Run human and mouse **separately**. For mouse downloads use
`data/file_mappings_borzoi_all.csv`, `--organism mouse`, and a separate output folder.
Use the raw or normalized matrix as `--truth` in the comparison command above.
Changing the gene/sample pool changes qnorm values; these scripts prepare new GT
and do not claim to recreate the historical matrices. Manifests record inputs.

## Notebooks

- [Prepare the support table from file mappings](notebooks/prepare_ontology_support_from_file_mappings.ipynb).
- [Explain ontology overlap and sample grouping](notebooks/ontology_mapping_explained.ipynb).
- [Model folds and gene TSS verification](notebooks/AG_versions.ipynb).

See [the official target workflow](../../docs/GENE_EXPRESSION_PAPER.md) for exon-coverage evaluation, and the [128 kb instructions](../../docs/RUN_128KB.md) for existing aglab inputs.
