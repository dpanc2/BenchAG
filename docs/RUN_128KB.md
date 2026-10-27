# AlphaGenome 128 kb expression run

Run on aglab in `/media/eternus1/nfs/projects/dpanchenko/BenchAG`.

- AlphaGenome client 0.7.0; model FOLD_1.
- Input length 131,072 bp (the client's legacy SEQUENCE_LENGTH_100KB constant).
- RNA_SEQ GeneMaskScorer, mean aggregation, 100,001 bp mask, gene-matched strand.
- Human and mouse, test then validation, sequential organism jobs.
- Human intervals: 2,781 test genes; 3,038 validation genes.
- Mouse intervals: 2,724 test genes; 2,998 validation genes.
- Inputs copied unmodified from anogena `/home/jovyan/dpanc/benchmarking/data`.
- Human truth file currently contains the historical 14-track panel. Predictions
  request all 285 mapped human tracks. Mouse truth contains 243 input columns;
  predictions request 68 mapped tracks. All predictions are retained regardless
  of overlap with truth. Broad-panel and non-qnorm scoring are separate analyses
  and must reuse these predictions with the correctly matched truth matrices.
- Original annotation release and exact truth preprocessing are unrecorded;
  configs state this explicitly. No new biological truth is constructed.
- Checkpoints saved every 5 genes. Genes absent from AlphaGenome annotation output are recorded in an outcomes
  JSON and skipped, allowing the next split to run. Genuine API/processing errors
  stop scoring after saving progress. Successful predictions resume from checkpoints.

Put the key in environment variable ALPHAGENOME_API_KEY or an ignored `.env` file
containing that assignment. BENCHAG_API_ENV may name another dotenv file. The
key is never included in configs or manifests. The driver removes proxy variables
for the AlphaGenome API child processes.

```bash
cd /media/eternus1/nfs/projects/dpanchenko/BenchAG
RUN_DIR="$PWD/results/runs/128kb_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$RUN_DIR"
nohup .venv/bin/python -u scripts/run_128kb.py --run-dir "$RUN_DIR" > "$RUN_DIR/driver.log" 2>&1 < /dev/null &
echo $! > "$RUN_DIR/driver.pid"
```

The driver creates run.json with PIDs, commands, logs and exit status. Human and
mouse logs are separate. Prediction manifests and checkpoints are under
results/expression/<run-id>/. A single-host lock prevents duplicate driver runs.
Before launch, validate both configs with `scripts/run_expression.py --check`.
Background launch is not completion; verify process state, logs and growing
checkpoint gene counts before reporting that inference has started.
