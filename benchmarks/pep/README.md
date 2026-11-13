# Enhancer–promoter interactions

Adapted from [NuriddinovMA/PEP-Benchmark](https://github.com/NuriddinovMA/PEP-Benchmark), commit `02ef6934e5cc5f7bf7870dea803f31b3298f94e4`.
The supplied CRISPR enhancer–gene pairs are included in `data/EP.GRCh38-c.tsv`.

Use Python 3.11, CUDA-enabled JAX, and access to the Hugging Face
`google/alphagenome-all-folds` checkpoint. Install `alphagenome==0.7.0`,
`alphagenome-research==0.3.0`, `pysam`, and `scikit-learn` in that environment.
Set up the [shared hg38 FASTA](../mpra/README.md#shared-genome), then run from
the repository root with a free GPU:

```bash
CUDA_VISIBLE_DEVICES=0 XLA_PYTHON_CLIENT_PREALLOCATE=false \
python benchmarks/pep/run.py --position-offset 1
```

Uses **ALL_FOLDS**, the shared `data/genomes/hg38.fa`, K562 total RNA-seq, and the recommended RNA-seq variant scorer. Each pair receives ten independently shuffled 2 kb enhancer replacements, with seed 42. Pair scores are the mean absolute raw effects across shuffles. Finished pairs are reused when restarting the same command.

Only pairs whose entire replacement fits in the promoter-centered 1 Mb input are scored: 6,606 of 10,356 pairs. Longer-distance pairs are excluded, rather than assigned zero. The final distance bin therefore covers only the portion that fits inside this input.

Results appear in `results/pep/all_folds/`: per-shuffle scores in `pairs/`, aggregated `predictions.tsv`, PR curves, and `metrics.json`. `--limit 2 --out results/pep/smoke` runs a small test.

Changes from the supplied scripts: local GPU inference; environment credentials rather than an embedded API key; indexed FASTA; seeded shuffles; resume support; working metrics for the draft’s distance groups. Labels come from `Regulated`, not `Significant`. `--position-offset 1` converts the zero-based FASTA start to AlphaGenome’s one-based variant position; `0` retains the original coordinate behavior. These corrections must be considered when comparing with the colleague’s saved values. auPRC uses raw scores. Pairs with missing gene or K562 scores are counted and excluded from auPRC.
