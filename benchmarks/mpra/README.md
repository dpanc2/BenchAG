# CAGI5 MPRA benchmark

Uses [tsoiegor/alphabench](https://github.com/tsoiegor/alphabench), pinned to
`58ed6e7a8d52cf88fdfea4125f5db21f723c5623`. The Enformer and Borzoi profiles
are ontology mappings; both run AlphaGenome ALL_FOLDS.

Settings: Tags > 10, source Position + 1, 1 Mb input, DNASE CenterMaskScorer
(width 501, DIFF_SUM), mean over matching tracks. Metrics are signed Pearson
and Spearman correlations with reporter activity, overall and per
(Element, Cell_Type). DNase variant effects serve as a proxy for MPRA effects.

## Setup

Use Python 3.11+, CUDA-enabled JAX, AlphaGenome Research, and access to ALL_FOLDS
weights. Clone the pinned package into an environment with those dependencies:

```bash
git clone https://github.com/tsoiegor/alphabench.git data/alphabench
cd data/alphabench
git checkout 58ed6e7a8d52cf88fdfea4125f5db21f723c5623
git apply ../../benchmarks/mpra/alphabench-splice-reference.patch
python -m pip install -e .
cd ../..
```

The patch supplies standard splice-site annotations needed by AlphaGenome
Research 0.3.0. Without them, variant inference fails at `jnp.asarray(None)`.
It does not change the DNase scorer or aggregation.

## Shared genome

The config uses `data/genomes/hg38.fa`, the same full hg38 reference as ATAC.
Place the reference there before running MPRA. On aglab it is already provided;
on a new machine you can download it from UCSC:

```bash
mkdir -p data/genomes
curl -L --fail --retry 3 \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  -o data/genomes/hg38.fa.gz
gzip -dc data/genomes/hg38.fa.gz > data/genomes/hg38.fa
```

## Run

From the repository root, download the original data, remove the single
reference-mismatched observation, and generate the run config:

```bash
mkdir -p data/mpra results/mpra
curl --fail --location 'https://zenodo.org/records/17453603/files/Kircher_GRCh38_ALL.tsv?download=1' \
  -o data/mpra/Kircher_GRCh38_ALL.tsv
python - <<'PYTHON'
from hashlib import md5, sha256
from pathlib import Path
import pandas as pd
import yaml

data = Path('data/mpra').resolve()
original = data / 'Kircher_GRCh38_ALL.tsv'
if md5(original.read_bytes()).hexdigest() != 'eacbb118660a3bd39575d88b1a03f7a7':
    raise ValueError('Unexpected CAGI5 input file')
source = pd.read_csv(original, sep='\t')
source.insert(0, 'source_row_id', source.index)
filtered = data / 'Kircher_GRCh38_REF_matched.tsv'
source.loc[source.source_row_id != 27570].to_csv(filtered, sep='\t', index=False)
config = yaml.safe_load(Path('benchmarks/mpra/cagi5.yaml').read_text())
config['data_dir'] = str(data)
config['output_dir'] = str(Path('results/mpra').resolve())
config['references']['hg38']['fasta_path'] = str(Path('data/genomes/hg38.fa').resolve())
settings = config['benchmarks']['cagi5_enformer_mapping']['config']
settings['variants_url'] = filtered.as_uri()
settings['variants_checksum'] = 'sha256:' + sha256(filtered.read_bytes()).hexdigest()
Path('results/mpra/cagi5.filtered.yaml').write_text(yaml.safe_dump(config))
PYTHON
alphabench download --config results/mpra/cagi5.filtered.yaml
alphabench process --config results/mpra/cagi5.filtered.yaml
CUDA_VISIBLE_DEVICES=0 alphabench run --config results/mpra/cagi5.filtered.yaml
```

Each mapping evaluates 29,034 observations across 22 assay groups.
Each run contains predictions, raw scores, and `metrics.json`. Use
`--resume <run-directory>` with the same config to reuse saved batches.

The run uses hg38 and the same 501 bp scorer for both tissue mappings.
The [AlphaGenome preprint methods](https://www.biorxiv.org/content/10.1101/2025.06.25.661532v1.full)
describe different scoring windows, a different reference assembly, and
specific locus exclusions for their comparison.
