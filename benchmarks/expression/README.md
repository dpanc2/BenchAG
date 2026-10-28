# AlphaGenome gene-expression benchmark

Start with your genes and original `file_mappings` CSV, match the samples to AlphaGenome,
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

Keep your original mapping table as `inputs/file_mappings.csv`. It needs an `id`
column with the sample identifiers used in the measured-expression table. Optional
`original_id` and `targets_identifier_base` columns provide alternate ENCODE
accessions. If `genome` or `assembly` is present, the script selects the requested
organism. Other columns, such as expression-file paths, can remain in the CSV.

```text
id,original_id,genome
ENCFF035CWS,ENCSR094GVZ,GRCh38
```

## 3. Find the ENCODE ontologies and build the support table

An ENCODE ID identifies an experimental file (`ENCFF…`) or experiment (`ENCSR…`).
An **ontology ID** identifies its tissue or cell type. AlphaGenome's output **name**
combines that ID and an assay, e.g. `UBERON:0006631 total RNA-seq`. Different assays
for the same ontology are different outputs.

Start from the original file mappings; no ontology table is needed:

```bash
python scripts/build_expression_support.py \
  --file-mappings inputs/file_mappings.csv --organism human
```

The script looks up each accession in ENCODE, follows its parent experiment when
needed, and retrieves the biosample ontology and assay. It then compares them with
AlphaGenome RNA-seq metadata and creates:

- `benchmarks/expression/data/qnorm_alphagenome_ontology_map.csv`: resolved ENCODE metadata.
- `benchmarks/expression/data/qnorm_alphagenome_support_summary.csv`: **one row per input sample**, with the AG comparison.
- `data/encode_metadata_cache.json`: raw ENCODE responses, reused on the next run.

AG metadata is reused from `alphagenome_human_rna_seq_metadata.csv` or exported
from the API if that file is absent. Use `--metadata` to choose another metadata
snapshot. The `qnorm` filenames follow the existing benchmark naming; this step
processes metadata only and does not normalize expression values.

Open the support-summary CSV and check:

| Column | What it tells you |
|---|---|
| `ontology_curie`, `biosample_name`, `assay` | Tissue/cell ontology and assay found in ENCODE |
| `alphagenome_supported` | Whether the ontology exists in AG, regardless of assay |
| `alphagenome_ontology_assay_supported` | **Whether the exact ontology AND assay name exists** |
| `alphagenome_ontology_name`, `alphagenome_strands` | Matched output name and available strands |
| `match_status`, `errors` | Unresolved lookups and unsupported names |

Use the **exact ontology–assay flag** to select prediction outputs. Unresolved
rows are retained, not silently dropped. Rerun the same command to reuse cached
ENCODE responses and fetch only missing ones.

Build the prediction/ground-truth JSON mappings from that resolved table:

```bash
python scripts/prepare_expression_mapping.py \
  --ontology-map benchmarks/expression/data/qnorm_alphagenome_ontology_map.csv \
  --metadata benchmarks/expression/data/alphagenome_human_rna_seq_metadata.csv \
  --organism human --out inputs/human_mapping
```

`track_to_first_id.json` gives one representative sample ID per AG name; use it
for prediction. `track_to_all_ids.json` keeps every matching sample ID; use it
when averaging measured samples for ground-truth comparison. Samples in the same
group share one prediction, but all their measurements can contribute to truth.
An ontology/assay match does not prove that the original experimental files are identical.

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

Repeat steps 2–5 with mouse genes (mm10), mouse file mappings, and
`--organism mouse` in the support-table, JSON-mapping, and prediction commands.
The generated ontology/support tables have `_mouse.csv` suffixes, and the AG
metadata file is `alphagenome_mouse_rna_seq_metadata.csv`. Use separate
mapping and output directories. No code edits are needed.

## Understand the matching and splits

Start with [prepare_ontology_support_from_file_mappings.ipynb](notebooks/prepare_ontology_support_from_file_mappings.ipynb)
to run the file-mapping → ENCODE → AG → support-table workflow cell by cell.
Then use [ontology_mapping_explained.ipynb](notebooks/ontology_mapping_explained.ipynb).
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
