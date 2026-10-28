# BenchAG

AlphaGenome benchmarks for gene expression, MPRA, and ATAC-seq.
The expression benchmark is ready; MPRA and ATAC-seq contributions are in progress.

## Start here

```bash
git clone https://github.com/dpanc2/BenchAG.git
cd BenchAG
bash scripts/setup.sh
source .venv/bin/activate
export ALPHAGENOME_API_KEY="YOUR_API_KEY"
```

Follow the [expression instructions](benchmarks/expression/README.md) to prepare
sample mappings, predict expression, and compare predictions with ground truth.
The sample-mapping CSVs are included in `data/`; large experimental files and API
keys are not included.

- [Prepare and inspect ontology mappings](benchmarks/expression/notebooks/prepare_ontology_support_from_file_mappings.ipynb)
- [MPRA contribution](benchmarks/mpra/README.md)
- [ATAC-seq contribution](benchmarks/atac/README.md)
- [How to contribute](CONTRIBUTING.md)
- [Data requirements](docs/DATA.md)

## Existing aglab runs

See [128 kb run instructions](docs/RUN_128KB.md) for the configured human and
mouse benchmark. Predictions and logs stay on the server under `results/`.

## Paper

The manuscript is being edited in the [shared Google Doc](https://docs.google.com/document/d/1ycHBmpq2U8bJkJv4r5r-CU4R524GsoIHuK6KXRdHNjQ/edit).
Source result records are retained under `paper/results/`.

## License

MIT; the imported expression code retains its upstream attribution.
See [LICENSE](LICENSE) and [source provenance](docs/UPSTREAM.md).
