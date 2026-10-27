# BenchAG: reproducible AlphaGenome benchmarks

A working repository for independent evaluation of AlphaGenome gene-expression,
MPRA, and ATAC-seq predictions. eQTL is optional and excluded from the initial
paper scope. Expression code is imported; other tasks await contributions.
No numerical benchmark results are asserted by this initial scaffold.

Primary workspace: `aglab:/media/eternus1/nfs/projects/dpanchenko/BenchAG`.

## Setup and one-command run

Python 3.11 is recommended for the pinned AlphaGenome client:

```bash
bash scripts/setup.sh
cp configs/expression.human.example.yaml configs/expression.human.yaml
# Edit REQUIRED metadata and data paths in the copied file.
export DATA_ROOT=/absolute/path/to/expression/data
export ALPHAGENOME_API_KEY=your_key
./run.sh expression --config configs/expression.human.yaml
```

Setup installs dependencies into `.venv`; `PYTHON_BOOTSTRAP` selects the setup
interpreter and `PYTHON_BIN` selects an existing runtime. Data and credentials are
external prerequisites. This repository does not download or invent ground truth.
For mouse use the mouse example; for 16 kb set sequence_length=16384 and
gene_mask_width=10001. Configuration never requires source-code edits.

Before spending API quota:

```bash
./run.sh expression --config configs/expression.human.yaml --check
```

The check needs pandas and PyYAML but no API key. It checks files, schema,
organism/strand labels, gene ID uniqueness, truth membership, mapping overlap,
and disjoint test/validation gene lists. It **does not certify** exclusion from
AlphaGenome or Borzoi training regions. Verify interval-level model fold overlap
separately before interpreting a result as held-out generalization.

## Inputs and outputs

See [data requirements](docs/DATA.md) for exact schemas and unresolved provenance.
Outputs are under `results/expression/<configuration-and-input-hash>/`:
resolved configuration, package/input manifest, per-split prediction CSVs,
checkpoints, and `Alphagenome_benchmark_summary.csv`. Repeat the same command to
resume a partial run. Failed genes stop final scoring so a partial table is not
silently presented as a completed benchmark. API calls can consume quota.

`mean_pearson_across_genes` averages per-track Pearson correlations over genes;
`mean_pearson_across_tracks` averages per-gene correlations over tracks. Legacy
aliases are retained. `deviation_r` is the imported selected-target pairwise
contrast metric with need_log=False; it is not a gene-mean-centered correlation.
Undefined correlations are omitted from the corresponding axis average. An
all-missing requested track remains visible in predictions. Missing selected
tracks make the historical specificity metric undefined; inspect coverage.

## Repository layout

- `benchmarks/expression/`: imported scripts, notebooks, and ontology mappings.
- `benchmarks/{mpra,atac,eqtl}/`: contribution slots with task-specific requirements.
- `configs/`: runnable configuration examples with required provenance fields.
- `scripts/`: common setup and expression entry point.
- `docs/`: data contract, source provenance, decisions, and contribution guidance.
- `paper/draft.md`: first scientific draft, with explicit missing-result markers.
- `paper/evidence.md`: claim status and evidence needed before submission.
- `results/`: generated local outputs, ignored by Git.

Run offline checks with `python3 -m unittest discover -s tests -v`.
Read [CONTRIBUTING.md](CONTRIBUTING.md) before adding another task.
Imported historical configs and Slurm scripts refer to the old GENA_LM layout;
they are retained for provenance, not recommended as launch commands.

The upstream MIT license is retained; see [provenance](docs/UPSTREAM.md).
A GitHub remote has not been created or selected yet.
