# BenchAG

AlphaGenome expression benchmark. The 128 kb setup on **aglab** is ready.

## Run

```bash
ssh aglab
cd /media/eternus1/nfs/projects/dpanchenko/BenchAG
RUN_DIR="results/runs/128kb_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$RUN_DIR"
nohup env BENCHAG_API_ENV="$PWD/.env.128kb" .venv/bin/python -u scripts/run_128kb.py --run-dir "$RUN_DIR" > "$RUN_DIR/driver.log" 2>&1 &
```

Runs human and mouse, test and validation. Missing-annotation genes are skipped;
existing checkpoints are reused. The API key stays in the ignored `.env.128kb` file.

## Results

- Predictions: `results/expression/<run-id>/predictions_{test,valid}_intervals.csv`.
- Logs and status: `$RUN_DIR/{human.log,mouse.log,run.json}`.
- Notebook gene/cell correlations: `results/correlations/` (qnorm and no-qnorm).

First-time setup: `bash scripts/setup.sh`, supply [input data](docs/DATA.md),
and edit the configs in `configs/`. Data and API keys are not included in Git.

[Full instructions](docs/EXPRESSION.md) · [128 kb details](docs/RUN_128KB.md) ·
[Comparison notebooks](benchmarks/expression/notebooks/) · [Contributing](CONTRIBUTING.md)
