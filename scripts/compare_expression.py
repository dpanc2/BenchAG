#!/usr/bin/env python3
"""Compare ontology-name predictions with measured expression, using the notebook metric."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

def pairwise_corr(x, y, axis):
    valid = np.isfinite(x) & np.isfinite(y)
    n = valid.sum(axis=axis)
    xx, yy = np.where(valid, x, 0.), np.where(valid, y, 0.)
    sx, sy = xx.sum(axis=axis), yy.sum(axis=axis)
    numerator = n * (xx * yy).sum(axis=axis) - sx * sy
    variance = (n * (xx * xx).sum(axis=axis) - sx ** 2) * (n * (yy * yy).sum(axis=axis) - sy ** 2)
    denominator = np.sqrt(np.maximum(variance, 0.))
    return np.divide(numerator, denominator, out=np.full(n.shape, np.nan), where=(n > 3) & (denominator > 0))

def group_values(frame, mapping, named_columns=False):
    values = {}
    for group, ids in mapping.items():
        if named_columns and group in frame:
            values[group] = pd.to_numeric(frame[group], errors='raise')
        else:
            present = [i for i in ids if i in frame]
            if present:
                values[group] = frame[present].astype(float).mean(axis=1)
    return pd.DataFrame(values, index=frame.index)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--predictions', required=True)
    p.add_argument('--truth', required=True, help='CSV: genes x samples (gene_id), or samples x genes (id)')
    p.add_argument('--mapping', required=True, help='track_to_all_ids.json')
    p.add_argument('--truth-transform', choices=['log2p1', 'none'], default='log2p1')
    p.add_argument('--out', required=True)
    args = p.parse_args()
    mapping = json.loads(Path(args.mapping).read_text())
    if not mapping or not all(isinstance(v, list) and v and all(isinstance(i, str) for i in v) for v in mapping.values()):
        p.error('Use the all-IDs mapping, not the representative-ID mapping')
    raw = pd.read_csv(args.truth)
    gene_columns = [c for c in raw if c.startswith(('ENSG', 'ENSMUSG'))]
    if gene_columns and ('id' in raw or 'gene_id' in raw):
        id_column = 'id' if 'id' in raw else 'gene_id'
        if raw[id_column].duplicated().any():p.error('Duplicate sample IDs in truth')
        truth = raw.set_index(id_column)[gene_columns].T.astype(float)
    elif 'gene_id' in raw:
        truth = raw.set_index('gene_id').astype(float)
    else:p.error('Truth needs gene_id rows, or sample id rows with Ensembl gene columns')
    prediction_raw = pd.read_csv(args.predictions)
    if 'gene_id' not in prediction_raw:p.error('Predictions need a gene_id column')
    pred = prediction_raw.set_index('gene_id')
    if truth.index.duplicated().any() or pred.index.duplicated().any():p.error('Duplicate gene IDs')
    if args.truth_transform == 'log2p1':
        if (truth < 0).any().any():p.error('log2p1 requires nonnegative truth values')
        truth = np.log2(truth + 1)
    # Transform individual experimental samples first, then average their group.
    truth = group_values(truth, mapping)
    pred = group_values(pred, mapping, named_columns=True)
    genes = truth.index.intersection(pred.index)
    groups = truth.columns.intersection(pred.columns)
    groups = [g for g in groups if np.isfinite(pred.loc[genes,g].to_numpy(dtype=float)).any()]
    if not len(genes) or not len(groups):p.error('No matching genes and available ontology-assay names')
    x, y = truth.loc[genes,groups].to_numpy(), pred.loc[genes,groups].to_numpy()
    per_track, per_gene = pairwise_corr(x,y,0), pairwise_corr(x,y,1)
    summary = {'gene_corr':float(np.nanmean(per_track)), 'cell_corr':float(np.nanmean(per_gene)),
               'n_genes':len(genes), 'n_ontology_assay_names':len(groups),
               'finite_gene_correlations':int(np.isfinite(per_gene).sum()),
               'finite_track_correlations':int(np.isfinite(per_track).sum()),
               'truth_transform':args.truth_transform}
    if not np.isfinite([summary['gene_corr'],summary['cell_corr']]).all():
        p.error('Correlations undefined: at least four finite, nonconstant observations are needed per axis')
    out = Path(args.out);out.mkdir(parents=True, exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame({'ontology_name':groups,'gene_corr':per_track}).to_csv(out/'per_ontology.csv',index=False)
    pd.DataFrame({'gene_id':genes,'cell_corr':per_gene}).to_csv(out/'per_gene.csv',index=False)
    (out/'inputs.json').write_text(json.dumps({str(Path(f).resolve()):hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in [args.predictions,args.truth,args.mapping]},indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
