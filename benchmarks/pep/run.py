"""Ten-shuffle enhancer ablation, adapted from NuriddinovMA/PEP-Benchmark."""
import argparse
import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd

p = argparse.ArgumentParser(__doc__)
p.add_argument('--genome', default='data/genomes/hg38.fa')
p.add_argument('--pairs', default=str(Path(__file__).parent / 'data/EP.GRCh38-c.tsv'))
p.add_argument('--out', default='results/pep/all_folds')
p.add_argument('--seed', type=int, default=42)
p.add_argument('--limit', type=int)
p.add_argument('--position-offset', type=int, choices=[0, 1], required=True)
p.add_argument('--prepare-only', action='store_true')
a = p.parse_args()
out = Path(a.out)
out.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(a.pairs, sep='\t')
data = data.drop_duplicates(['measuredGeneEnsemblId', 'startTSS', 'chromStart'], keep='last')
eligible = data[(data.startTSS - 524288 < data.center) &
                (data.startTSS + 524288 > data.center + 2000)].copy()
eligible.to_csv(out / 'eligible_pairs.tsv', sep='\t', index=False)
(out / 'settings.json').write_text(json.dumps({**vars(a), 'checkpoint': 'all_folds',
    'input_pairs': len(data), 'eligible_pairs': len(eligible), 'shuffles': 10,
    'label': 'Regulated', 'ontology': 'EFO:0002067', 'assay': 'total RNA-seq'}, indent=2))
print(f'{len(eligible)}/{len(data)} pairs fit the 1 Mb input; ten shuffles each', flush=True)
if a.prepare_only:
    raise SystemExit

import pysam
from alphagenome.data import genome
from alphagenome.models import variant_scorers
from alphagenome_research.model import dna_model

if not Path(a.genome + '.fai').exists():
    pysam.faidx(a.genome)

settings = dict(dna_model.default_organism_settings())
human = dna_model.Organism.HOMO_SAPIENS
settings[human] = dataclasses.replace(settings[human], fasta_path=a.genome,
                                     calibration_path=None, pas_feather_path=None)
settings[dna_model.Organism.MUS_MUSCULUS] = dna_model.OrganismSettings()
model = dna_model.create_from_huggingface('all_folds', organism_settings=settings)
fasta = pysam.FastaFile(a.genome)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
cache = out / 'pairs'
cache.mkdir(exist_ok=True)
rows = list(eligible.sort_values(['measuredGeneEnsemblId', 'startTSS', 'chromStart']).itertuples())
if a.limit:
    rows = rows[:a.limit]
for i, r in enumerate(rows):
    dest = cache / f'{i:05d}.tsv'
    if dest.exists():
        continue
    reference = fasta.fetch(r.chrom, r.center, r.center + 2000).upper()
    rng = np.random.default_rng(np.random.SeedSequence([a.seed, i]))
    interval = genome.Interval(r.chrTSS, r.startTSS - 524288, r.startTSS + 524288)
    scores = []
    for j in range(10):
        variant = genome.Variant(chromosome=r.chrom, position=r.center + a.position_offset,
            reference_bases=reference, alternate_bases=''.join(rng.permutation(list(reference))))
        result = model.score_variants([interval], [variant], variant_scorers=[scorer],
                                      max_workers=1, progress_bar=False)
        df = variant_scorers.tidy_scores(result)
        df = df[(df.gene_id == r.measuredGeneEnsemblId) &
                (df.ontology_curie == 'EFO:0002067') & (df['Assay title'] == 'total RNA-seq')]
        if df.empty:
            print(f'{i + 1}/{len(rows)} {r.measuredGeneEnsemblId}: missing gene/K562 score', flush=True)
            scores = []
            break
        scores.append({'shuffle': j, 'raw_score': float(df.raw_score.mean())})
    pd.DataFrame(scores, columns=['shuffle', 'raw_score']).to_csv(dest.with_suffix('.tmp'), sep='\t', index=False)
    dest.with_suffix('.tmp').replace(dest)
    print(f'{i + 1}/{len(rows)} {r.measuredGeneEnsemblId}', flush=True)

from sklearn.metrics import auc, precision_recall_curve
predictions = []
for i, r in enumerate(rows):
    df = pd.read_csv(cache / f'{i:05d}.tsv', sep='\t')
    predictions.append({'gene_id': r.measuredGeneEnsemblId, 'chrom': r.chrom,
        'center': r.center, 'distance': r.distanceToTSS,
        'label': str(r.Regulated).upper() == 'TRUE',
        'score': float(df.raw_score.abs().mean()) if len(df) == 10 else np.nan})
predictions = pd.DataFrame(predictions)
predictions.to_csv(out / 'predictions.tsv', sep='\t', index=False)
metrics = []
for name, mask in [('all', predictions.distance >= 0),
                   ('0-10kb', predictions.distance < 10000),
                   ('10-100kb', predictions.distance.between(10000, 100000, inclusive='left')),
                   ('100-2500kb', predictions.distance.between(100000, 2500000, inclusive='left'))]:
    candidates = predictions[mask]
    subset = candidates.dropna(subset=['score'])
    value = None
    if subset.label.nunique() == 2:
        precision, recall, _ = precision_recall_curve(subset.label, subset.score)
        value = float(auc(recall, precision))
        pd.DataFrame({'precision': precision, 'recall': recall}).to_csv(out / f'pr_{name}.csv', index=False)
    metrics.append({'distance': name, 'n_pairs': len(subset), 'n_missing': len(candidates) - len(subset), 'n_positive': int(subset.label.sum()), 'auprc': value})
(out / 'metrics.json').write_text(json.dumps(metrics, indent=2))
print(json.dumps(metrics, indent=2), flush=True)
