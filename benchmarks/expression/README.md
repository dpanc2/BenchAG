# AlphaGenome gene-expression benchmark

Start with your genes and ENCODE sample IDs, match the samples to AlphaGenome,
then predict expression and compare it with measured values. Human comes first
below; the same commands support mouse.

**Run every command from the BenchAG repository root.**

## 1. Set up BenchAG

Replace the clone URL below with the URL of your BenchAG repository.

```bash
git clone <your-BenchAG-repository-URL> BenchAG
cd BenchAG
conda create -n benchag python=3.11 pip -y
conda activate benchag
python -m pip install -r requirements.txt
export ALPHAGENOME_API_KEY="YOUR_API_KEY"
```

Already have BenchAG and its environment? Activate that environment and start at
step 2. Keep your API key outside Git. The requirements install the tested client
and scientific packages together; record the installed client version when
reporting results (`python -m pip show alphagenome`). Strand selection is handled
by the prediction command, so you do not edit the Python source for each strand.

## 2. Prepare genes and samples

Create an `inputs/` folder for your files.

Create `inputs/genes.tsv` with these columns (tab-separated):

```text
gene_id             chrom   TSS       gene_strand
ENSG00000001036.13   chr6    1000000   +
```

The row illustrates the format: replace it with your gene's real chromosome,
**zero-based TSS**, and strand (`+` or `-`). Use one TSS per gene and the matching
genome assembly: **hg38 for human** or **mm10 for mouse**. Gene IDs may include a
version suffix. No reference FASTA or ground-truth file is needed for prediction.

Create `inputs/samples.txt` with one ENCODE file (`ENCFF…`) or experiment
(`ENCSR…`) ID per line. These IDs should be the same identifiers used in your
measured-expression table. A CSV/TSV with an `id` column also works.

## 3. Match samples to AlphaGenome output names

An ENCODE ID identifies an experimental file or experiment. An **ontology ID**
identifies its tissue or cell type, such as `CL:0000047`. AlphaGenome's output
**name** combines an ontology ID and an assay, for example
`UBERON:0006631 total RNA-seq`. We call this an **ontology–assay name**: the same
ontology with a different RNA-seq assay is a different output.

```bash
python scripts/prepare_expression_mapping.py \
  --samples inputs/samples.txt --organism human --out inputs/human_mapping
```

This command gets the tissue/cell ontology and assay from ENCODE, exports
AlphaGenome RNA-seq metadata, and matches both fields. Inspect
`inputs/human_mapping/sample_matching.csv` before prediction: it lists matched,
unresolved, and unsupported samples. Unmatched samples are not guessed.

It also creates two JSON files:

- `track_to_first_id.json`: each AlphaGenome name → one representative sample ID;
  pass this file to the prediction command.
- `track_to_all_ids.json`: each name → all matching sample IDs; use this file to
  average experimental samples when comparing with ground truth.

Samples in the same group share one predicted value; the representative ID does
not mean that the other measurements are discarded. Matching an ontology and
assay does **not** mean the model and benchmark used identical experimental files.

## 4. Run human predictions

```bash
python scripts/predict_expression.py \
  --genes inputs/genes.tsv \
  --mapping inputs/human_mapping/track_to_first_id.json \
  --organism human --window 128kb --out results/my_human_run
```

Choose `--window 16kb`, `128kb`, or `1mb`. The command selects the organism,
sequence length, scoring width, and gene strand automatically. The model defaults
to `FOLD_1`; use `--model-version` to select another version deliberately.
Add `--check` to validate the input files without making prediction requests.

The output directory contains:

- `predictions.csv`: **one row per gene, one column per ontology–assay name**.
- `checkpoint.csv` and `checkpoint.csv.outcomes.json`: saved progress and skipped
  genes. Rerun the unchanged command to resume.
- `prediction_manifest.json`: input hashes, settings, and client version.

Genes absent from the annotation output are skipped and recorded. An all-`NaN`
column means the requested output was unavailable under the runner's strand
selection; it is not a zero-expression prediction. Genuine request errors are
recorded separately. Do not assume that every requested name has usable values.

For a background run, put `nohup` before the command and append
`> results/my_human_run.log 2>&1 < /dev/null &`.

## 5. Compare with ground truth

Use a measured-expression CSV with `gene_id` rows and sample-ID columns. The
transposed layout (`id` rows and Ensembl gene columns) is also accepted.

```bash
python scripts/compare_expression.py \
  --predictions results/my_human_run/predictions.csv \
  --truth inputs/truth.csv \
  --mapping inputs/human_mapping/track_to_all_ids.json \
  --out results/my_human_run/comparison
```

By default, the command applies `log2(x + 1)` to each experimental sample, then
averages samples with the same ontology–assay name. If the truth is already on the
scale you want to evaluate, add `--truth-transform none`. To compare qnorm and
non-qnorm targets, supply the corresponding truth files separately; this command
does not create quantile-normalized data.

`comparison/summary.json` reports:

- **gene correlation**: Pearson across genes, averaged over output names;
- **cell correlation**: Pearson across output names, averaged over genes.

Only matching genes/groups and finite paired values are evaluated. Constant or
insufficient vectors are excluded, and the evaluated counts are saved.

## Mouse

Repeat steps 2–5 with mouse genes (mm10), mouse sample IDs, and
`--organism mouse` in **both** the mapping and prediction commands. Use separate
mapping and output directories. No code edits are needed.

## Understand the matching and splits

Start with [ontology_mapping_explained.ipynb](notebooks/ontology_mapping_explained.ipynb).
It recomputes the human and mouse overlap counts and explains tracks, strands,
ontology IDs, ontology–assay names, and the two JSON mappings. The human reference
benchmark has **815 RNA-seq tracks**, covering tissues and cell lines; these are
not 815 distinct cell types. Its 770 matched tracks form 285 shared ontology–assay
names with 240 distinct ontology IDs.

[AG_versions.ipynb](notebooks/AG_versions.ipynb) shows model folds. Borzoi test
(fold3) is AlphaGenome FOLD_1 validation; Borzoi validation (fold4) is AlphaGenome
test. The existing benchmark tables keep Borzoi's split names.

The [human](notebooks/gena_ag_gt_ontology_compare.ipynb) and
[mouse](notebooks/gena_ag_gt_ontology_compare_mm10.ipynb) comparison notebooks
retain the historical benchmark analysis. For your own samples and the
ontology-name output above, use step 5.
