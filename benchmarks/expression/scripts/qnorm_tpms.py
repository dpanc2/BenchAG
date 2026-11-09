#!/usr/bin/env python3
"""Quantile-normalize a gene-by-sample TPM CSV for one species."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
import qnorm


def normalize(tpms, out):
    frame = pd.read_csv(tpms, dtype={'gene_id':str})
    if 'gene_id' not in frame:
        raise ValueError('Input needs a gene_id column and numeric sample columns')
    frame = frame.set_index('gene_id')
    if frame.empty or frame.index.has_duplicates or frame.index.isna().any():
        raise ValueError('Input needs unique, nonempty gene IDs')
    # Accept the split label present in older GT matrices, without normalizing it.
    splits = frame.pop('split') if 'split' in frame else None
    frame = frame.apply(pd.to_numeric, errors='raise').astype(float)
    if frame.shape[1] == 0 or not np.isfinite(frame.to_numpy()).all() or (frame < 0).any().any():
        raise ValueError('TPMs must be finite, nonnegative numbers without missing values')
    result = qnorm.quantile_normalize(frame, axis=1, ncpus=1)
    if splits is not None:
        result['split'] = splits
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index_label='gene_id')
    manifest = {'input':str(tpms), 'input_sha256':hashlib.sha256(Path(tpms).read_bytes()).hexdigest(),
                'genes':len(frame), 'samples':frame.shape[1],
                'normalization':'raw TPM; qnorm.quantile_normalize(axis=1); no log or filtering',
                'versions':{p:importlib.metadata.version(p) for p in ['qnorm','numpy','pandas']}}
    out.with_suffix('.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Saved {len(frame)} genes, {frame.shape[1]} samples to {out}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tpms', required=True, help='One species; gene_id column followed by samples')
    parser.add_argument('--out', required=True)
    normalize(**vars(parser.parse_args()))
