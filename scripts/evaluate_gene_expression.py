#!/usr/bin/env python3
"""Three gene-expression correlations from AlphaGenome Fig. 2d (2025 preprint).

Accepts gene x track matrices, optionally grouping TPM samples by ontology/assay.
This implements the reported metric, not a claim of identical evaluation data.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def quantile_normalize(x):
    """Equalize column distributions; ties receive the mean of tied rank values.

    Predicted and observed matrices are normalized independently. Missing values
    must be resolved by the caller; they must never become zero expression.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or not np.isfinite(x).all() or not x.size:
        raise ValueError('Quantile normalization needs a nonempty finite matrix')
    target = np.sort(x, axis=0).mean(axis=1)
    out = np.empty_like(x)
    for j in range(x.shape[1]):
        low = rankdata(x[:, j], method='min').astype(int) - 1
        high = rankdata(x[:, j], method='max').astype(int)
        prefix = np.r_[0., np.cumsum(target)]
        out[:, j] = (prefix[high] - prefix[low]) / (high - low)
    return out


def normalize_expression(x):
    q = quantile_normalize(x)
    return q - q.mean(axis=1, keepdims=True)


def correlations(x, y, axis, min_pairs=3):
    valid = np.isfinite(x) & np.isfinite(y)
    n = valid.sum(axis=axis)
    mx = np.divide(np.where(valid, x, 0).sum(axis=axis), n,
                   out=np.zeros(n.shape), where=n > 0)
    my = np.divide(np.where(valid, y, 0).sum(axis=axis), n,
                   out=np.zeros(n.shape), where=n > 0)
    dx = np.where(valid, x - np.expand_dims(mx, axis), 0.)
    dy = np.where(valid, y - np.expand_dims(my, axis), 0.)
    denominator = np.sqrt((dx * dx).sum(axis=axis) * (dy * dy).sum(axis=axis))
    r = np.divide((dx * dy).sum(axis=axis), denominator,
                  out=np.full(n.shape, np.nan), where=(n >= min_pairs) & (denominator > 0))
    return np.clip(r, -1, 1), n


def evaluate(prediction, observed):
    """Evaluate the same complete gene/track rectangle for all three metrics."""
    if not prediction.index.is_unique or not observed.index.is_unique:
        raise ValueError('Gene IDs must be unique')
    if not prediction.columns.is_unique or not observed.columns.is_unique:
        raise ValueError('Track IDs must be unique')
    genes = observed.index.intersection(prediction.index, sort=False)
    tracks = observed.columns.intersection(prediction.columns, sort=False)
    p = prediction.loc[genes, tracks].to_numpy(dtype=float)
    o = observed.loc[genes, tracks].to_numpy(dtype=float)
    available = (np.isfinite(p) & np.isfinite(o)).any(axis=0)
    dropped_tracks = tracks[~available].tolist()
    tracks, p, o = tracks[available], p[:, available], o[:, available]
    complete = (np.isfinite(p) & np.isfinite(o)).all(axis=1)
    dropped_genes = genes[~complete].tolist()
    genes, p, o = genes[complete], p[complete], o[complete]
    if len(genes) < 3 or len(tracks) < 3:
        raise ValueError('Need at least three complete genes and three available tracks')
    pn, on = normalize_expression(p), normalize_expression(o)
    rows = []
    for metric, a, b, axis, ids in [
        ('cross_gene_raw', p, o, 0, tracks),
        ('cross_gene_normalized', pn, on, 0, tracks),
        ('cross_track_normalized', pn, on, 1, genes),
    ]:
        rs, ns = correlations(a, b, axis)
        rows.extend({'metric': metric, 'id': str(key), 'pearson_r': r, 'n_pairs': int(n)}
                    for key, r, n in zip(ids, rs, ns))
    scores = pd.DataFrame(rows)
    summary = scores.groupby('metric', sort=False).pearson_r.agg(['mean', 'median', 'count']).reset_index()
    summary['n_genes'] = len(genes)
    summary['n_tracks'] = len(tracks)
    if summary['count'].eq(0).any():
        raise ValueError('A metric is undefined for every gene/track (constant inputs)')
    coverage = {'aligned_genes': len(complete), 'aligned_tracks': len(available),
                'dropped_all_missing_tracks': dropped_tracks,
                'dropped_incomplete_genes': dropped_genes,
                'used_genes': genes.tolist(), 'used_tracks': tracks.tolist(),
                'prediction_only_genes': prediction.index.difference(observed.index).tolist(),
                'observed_only_genes': observed.index.difference(prediction.index).tolist(),
                'prediction_only_tracks': prediction.columns.difference(observed.columns).tolist(),
                'observed_only_tracks': observed.columns.difference(prediction.columns).tolist()}
    return summary, scores, coverage


def read_matrix(path):
    with open(path, newline='') as handle:
        header = next(csv.reader(handle), [])
    if len(header) != len(set(header)):
        raise ValueError(f'{path}: duplicate column names')
    frame = pd.read_csv(path, index_col=0, dtype={'gene_id': str}, low_memory=False)
    if frame.index.name != 'gene_id':
        raise ValueError(f'{path}: expected gene_id as the first column')
    if frame.empty or frame.index.hasnans or not frame.index.is_unique:
        raise ValueError(f'{path}: empty matrix, missing or duplicate gene IDs')
    return frame.astype(float)


def transform(frame, name):
    if name == 'none':
        return frame
    if (frame < 0).any().any():
        raise ValueError('Log transform requires nonnegative input')
    return np.log1p(frame)


def group_tracks(frame, mapping, named=False):
    # Require every replicate; a partially missing group is not silently averaged.
    groups = {}
    for name, samples in mapping.items():
        if named and name in frame:
            groups[name] = frame[name]
        elif all(s in frame for s in samples):
            groups[name] = frame[samples].mean(axis=1, skipna=False)
    return pd.DataFrame(groups, index=frame.index)


def plot_results(scores, destination, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    metrics = ['cross_gene_raw', 'cross_gene_normalized', 'cross_track_normalized']
    fig, ax = plt.subplots(figsize=(6.3, 4.4))
    for pos, metric in enumerate(metrics, 1):
        values = scores.loc[scores.metric == metric, 'pearson_r'].dropna().to_numpy()
        if len(values) > 1 and np.ptp(values) > 1e-12:
            violin = ax.violinplot(values, positions=[pos], widths=.65, showextrema=False)
            violin['bodies'][0].set_facecolor('#2a9d8f')
            violin['bodies'][0].set_alpha(.65)
        ax.scatter(pos, values.mean(), color='#e76f51', zorder=3)
        ax.annotate(f'{values.mean():.3f}', (pos, values.mean()), xytext=(8, 8), textcoords='offset points')
    ax.set_xticks([1, 2, 3], ['Raw\nCross-gene', 'Normalized\nCross-gene', 'Normalized\nCross-track'])
    ax.set_ylabel('Pearson r per track or gene')
    ax.set_ylim(-1.02, 1.08)
    ax.set_title(title, fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions', required=True)
    parser.add_argument('--observed', required=True)
    parser.add_argument('--mapping', help='Optional ontology-assay -> sample IDs JSON for TPM comparison')
    parser.add_argument('--prediction-transform', choices=['log1p', 'none'], required=True)
    parser.add_argument('--observed-transform', choices=['log1p', 'none'], required=True)
    parser.add_argument('--input-kind', choices=['exon_coverage', 'tpm_proxy', 'other'], required=True)
    parser.add_argument('--label', default='Gene-level RNA-seq prediction')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    pred = transform(read_matrix(args.predictions), args.prediction_transform)
    obs = transform(read_matrix(args.observed), args.observed_transform)
    if args.mapping:
        mapping = json.loads(Path(args.mapping).read_text())
        if not mapping or not all(isinstance(v, list) and v and all(isinstance(s, str) for s in v)
                                  for v in mapping.values()):
            raise ValueError('Mapping must contain nonempty lists of sample IDs')
        # Experimental replicates are transformed BEFORE averaging, like BenchAG.
        obs = group_tracks(obs, mapping)
        pred = group_tracks(pred, {k: [v[0]] for k, v in mapping.items()}, named=True)
    summary, scores, coverage = evaluate(pred, obs)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'manifest.json').exists():
        raise ValueError('Output already contains an evaluation; use a new --out directory')
    summary.to_csv(out / 'summary.csv', index=False)
    scores.to_csv(out / 'correlations.csv', index=False)
    (out / 'coverage.json').write_text(json.dumps(coverage, indent=2) + '\n')
    inputs = [args.predictions, args.observed] + ([args.mapping] if args.mapping else [])
    manifest = {'arguments': vars(args), 'normalization': 'independent quantile normalization per matrix; gene mean centering',
                'ties': 'mean target value over tied ranks', 'missing': 'drop all-missing tracks, then incomplete genes',
                'paper_reproduction_certified': False,
                'source': 'https://storage.googleapis.com/deepmind-media/papers/alphagenome.pdf',
                'sha256': {str(Path(f).resolve()): hashlib.sha256(Path(f).read_bytes()).hexdigest()
                           for f in inputs + [__file__]}}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    plot_results(scores, out / 'gene_expression.png', args.label + '\n' + args.input_kind)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
