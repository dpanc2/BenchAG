# IMR-90 ATAC-seq benchmark

Reproduces [GalinaK413/atac_bench_example](https://github.com/GalinaK413/atac_bench_example),
commit `9e6a91157702b80efbc3b43f23c06e21f7bae856`. The wrappers call the colleague's
prediction and metric functions, adding explicit model versions and saved profiles
for resuming. Results stay under `results/atac/`.

## Inputs

Clone the pinned source into `data/atac/source`. Download these ENCODE files into
`data/atac/` (URLs: `https://www.encodeproject.org/files/<accession>/@@download/<filename>`):

- `ENCFF650PHA.bigWig`: observed IMR-90 ATAC signal.
- `ENCFF113GSV.tar.gz`: ChromBPNet models; extract the fold-0 no-bias and scaled-bias H5 files.
- `ENCFF802KJC.tar.gz`: original training/test regions, used to verify the supplied BED.

Their checked MD5 values are recorded in `results/atac/input_checksums.json`.
Use the source's `alpha/fold1_test_regions.bed` for every model: **13,812 peak
records, all from ChromBPNet fold-0 test peaks**. Duplicate evaluation windows
are retained in metrics; API predictions are cached by unique window.

## AlphaGenome

Use the repository environment plus `python -m pip install pyBigWig` and set
`ALPHAGENOME_API_KEY`. Run from the BenchAG root:

```bash
python benchmarks/atac/run_alphagenome.py --model-version FOLD_1 --out results/atac/fold_1
python benchmarks/atac/run_alphagenome.py --model-version FOLD_0 --out results/atac/fold_0
```

Use `--limit 2` with a separate output directory for a smoke test. Rerun the same
command to reuse saved region arrays. `--workers 2` controls concurrent windows.
Settings match the source: 1 Mb input, IMR-90 ontology EFO:0001196, ATAC output,
1 kb evaluation windows centered at peak midpoints, forward/reverse averaging,
1 bp global Pearson, log1p total-count Pearson, and mean Jensen–Shannon distance.
The unused upstream torch import is omitted; no PyTorch installation is needed.

**Fold audit:** 13,550 of the 13,812 evaluation windows overlap FOLD_0 training
regions. None overlap FOLD_1 training regions. Therefore FOLD_0 reproduces the
script as a diagnostic; FOLD_1 is the held-out comparison recommended by the
source README. This distinction is recorded with the results.

## ChromBPNet

Use a separate CPU environment with TensorFlow CPU 2.15.1, TensorFlow Probability
0.23.0, numpy, pandas, scipy, matplotlib, tqdm, pyBigWig and pysam. Install
ChromBPNet 0.1.7 with `pip install --no-deps chrombpnet==0.1.7` to retain these
TensorFlow versions. Combine the two fold-0 models using the
[cited reformat helper](https://github.com/kundajelab/chrombpnet/blob/master/chrombpnet/helpers/postprocessing/reformat_chrombpnet_h5.py),
writing `data/atac/models/chrombpnet_recompiled.h5`.

```bash
data/atac/chrombpnet_env/bin/python benchmarks/atac/run_chrombpnet.py \
  --genome data/atac/hg38_chr1_chr3_chr6.fa --out results/atac/chrombpnet
```

The local reference contains the unchanged hg38 sequences for chr1, chr3 and
chr6, the only chromosomes in the supplied BED. Its origin and per-chromosome
sequence hashes are saved in `data/atac/reference_subset.json`. A full hg38 FASTA
can also be passed. The wrapper indexes the FASTA if needed.

Predictions use summit-centered 2,114 bp inputs and 1,000 bp outputs, with
reverse-complement averaging and overlap averaging into a BigWig. Evaluation
uses midpoint-centered windows, as in the source; uncovered prediction bases
are treated as zero by its metrics. The wrapper reports their number separately.

## Saved results

Each model output directory contains `metrics.json`, `summary.csv`, an input
manifest, and cached region profiles. ChromBPNet also saves its predicted BigWig.
AlphaGenome API runs and ChromBPNet CPU inference run on aglab. The source
paper values are comparison targets, not results of these new runs.
