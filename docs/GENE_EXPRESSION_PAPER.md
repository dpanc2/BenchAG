# Gene expression against released AlphaGenome targets

Evaluate FOLD_1 predictions against the released human or mouse RNA-seq
TFRecords. The workflow compares mean exon coverage on the same bases in
prediction and target, selects tracks by gene strand, and preserves missing
targets. It reports gene, track, and gene-by-context correlations. Exact
agreement with the original figure also depends on the authors' evaluator
and fold aggregation.

Run from the repository root on Linux with Python 3.11 and an AlphaGenome
API key. The pilot reads three gene-bearing windows. Full runs stream tens
of gigabytes and make multiple API calls per window.

```bash
bash scripts/setup.sh
source .venv/bin/activate
python -m pip install tensorflow-cpu==2.20.0 google-crc32c
export ALPHAGENOME_API_KEY='YOUR_API_KEY'

mkdir -p data/official_targets
curl --fail --location https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz \
  -o data/official_targets/gencode.v46.primary_assembly.annotation.gtf.gz
python scripts/run_official_expression.py \
  --gcs-prefix v1/train/FOLD_1/HOMO_SAPIENS/TEST/RNA_SEQ/ \
  --metadata benchmarks/expression/data/target_human_metadata.csv \
  --gtf data/official_targets/gencode.v46.primary_assembly.annotation.gtf.gz \
  --organism HOMO_SAPIENS --model-version FOLD_1 --api-env .env \
  --max-records 3 --ontology-batch-size 64 \
  --out results/official_expression/human_pilot
```

For a full run, set `--max-records 0` and use a new output directory.
For mouse, replace `HOMO_SAPIENS` with `MUS_MUSCULUS` in both places, use
`target_mouse_metadata.csv`, and download
[GENCODE M23](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz)
for mm10. Use `--records` instead of `--gcs-prefix` for local TFRecords.
The metadata CSVs preserve the released channel order, including padding.
The [official dataset loader](https://github.com/google-deepmind/alphagenome_research/blob/main/src/alphagenome_research/io/dataset.py)
describes the source records.

Rerun the same command to resume. The output contains paired exon-coverage
matrices, metrics, and plots. To evaluate previously generated matrices, use
`python scripts/evaluate_gene_expression.py --help`. To construct matrices from
BigWig tracks, install `pyBigWig` and use
`python scripts/aggregate_gene_exons.py --help`.
