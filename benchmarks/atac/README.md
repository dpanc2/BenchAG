# ATAC-seq benchmark

Compare AlphaGenome FOLD_0, FOLD_1 and ChromBPNet on IMR-90 ATAC-seq data.
Requires Linux, Python 3.10/3.11 and an AlphaGenome API key. No GPU is needed.

## 1. Open the project and load your key

From the BenchAG repository root:

```bash
bash scripts/setup.sh
source .venv/bin/activate
export ALPHAGENOME_API_KEY='paste-your-key-here'
```

## 2. Prepare the BED files once

Open [`prepare_inputs.ipynb`](prepare_inputs.ipynb), select the BenchAG `.venv`
kernel and run all cells. It creates both BED files in `data/atac/`.
If needed, install its dependencies first:

```bash
pip install alphagenome==0.7.0 pyranges==0.1.4 pandas
```

## 3. Try two regions first

```bash
bash benchmarks/atac/run_benchmark.sh results/atac/quick_test 2
```

This installs the benchmark dependencies, downloads missing signal/model files
and tests all three models on two regions. Prepare the BEDs in the notebook first.
A small results table appears when it finishes. These two-region metrics are
only a test, not the final benchmark results.

## 4. Run the full benchmark

```bash
mkdir -p results/atac
nohup bash benchmarks/atac/run_benchmark.sh results/atac/full_run \
  > results/atac/full_run.log 2>&1 &
```

The job continues if you disconnect. Follow setup progress with:

```bash
tail -f results/atac/full_run.log
```

Once prediction starts, each model writes its own progress log:

```bash
tail -f results/atac/full_run/fold_1/run.log
```

The other logs are `fold_0/run.log` and `chrombpnet/run.log` in the same folder.
Press **Ctrl+C** to stop watching a log; the background job continues.

## Results and restarting

The final table is `results/atac/full_run/summary.csv`. Each model folder
contains `metrics.json`, `summary.csv` and saved prediction profiles.
The metrics are global Pearson correlation, log-count Pearson correlation,
and mean Jensen–Shannon distance.

If a run fails, read its log and rerun the same command to reuse saved profiles.
For fresh predictions, change `full_run` to a new folder name. Do not reuse the
`quick_test` folder for a full run.

## Data and BED files

Signal and model files are downloaded into `data/atac/`:

- [ENCFF650PHA](https://www.encodeproject.org/files/ENCFF650PHA/): observed ATAC signal.
- [ENCFF113GSV](https://www.encodeproject.org/files/ENCFF113GSV/): ChromBPNet model weights.
- [ENCFF802KJC](https://www.encodeproject.org/files/ENCFF802KJC/): original IMR-90 peak regions.

The full human reference is **`data/genomes/hg38.fa`**, shared with MPRA.
ATAC reuses it if present; otherwise it downloads the full UCSC hg38 genome.
There is no separate ATAC chromosome-subset reference.

[`prepare_inputs.ipynb`](prepare_inputs.ipynb) follows
[this Colab](https://colab.research.google.com/drive/1qjR2W1LpRbrb0zmGescaiUGe7Ft1k0kN):
it removes peaks overlapping AG training regions expanded to 1 Mb, and saves
`fold0_test_regions.bed` (11,860 peaks) and `fold1_test_regions.bed` (13,812 peaks).
The Colab's original archive, ENCFF574LOE, is for a different experiment.

Our comparison uses **the same 13,812 FOLD_1-filtered peaks for all three models**.
FOLD_1 is held out; FOLD_0 is a diagnostic because most of these regions overlap
its training partitions. Prediction and metric code is included in `alpha/` and
`chrombpnet/`, adapted from [the colleague's benchmark](https://github.com/GalinaK413/atac_bench_example/tree/9e6a91157702b80efbc3b43f23c06e21f7bae856).
