#!/usr/bin/env python3
"""Prepare ENCODE TPM and qnorm targets; adapted from minjaf/ExpressionPredict.

Source: 6c10d3b7f2022488f40c6d3c52ddca338756cea5,
src/download_borzoi.py, src/qnorm_tpm_data.py and src/utils.py.
Uses recorded quantification accessions, so no BAM downloads are needed.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests


def prepare(file_mappings, organism, annotations, intervals_prefix, out):
    import qnorm

    out = Path(out)
    cache = out / 'encode'
    cache.mkdir(parents=True, exist_ok=True)
    assemblies = {'human': ['GRCh38', 'hg38'], 'mouse': ['mm10']}[organism]
    mappings = pd.concat([pd.read_csv(p, dtype=str) for p in file_mappings], ignore_index=True)
    mappings = mappings[mappings.genome.isin(assemblies)].copy()
    if mappings.empty:
        raise ValueError(f'No {organism} samples in file mappings')
    if mappings.encode_tsv.isna().any():
        raise ValueError('Every sample must have an encode_tsv quantification accession')
    for sample, group in mappings.groupby('id'):
        if group.encode_tsv.nunique() != 1:
            raise ValueError(f'Conflicting quantification accessions for {sample}')
    mappings = mappings.drop_duplicates('id')
    series = []
    sources = []
    for row in mappings.itertuples():
        accession = str(row.encode_tsv).removesuffix('.tsv')
        if not accession.startswith('ENCFF') or not accession[5:].isalnum():
            raise ValueError(f'Invalid ENCODE accession: {accession}')
        url = f'https://www.encodeproject.org/files/{accession}/@@download/{accession}.tsv'
        path = cache / f'{accession}.tsv'
        if not path.exists():
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            temporary = path.with_suffix('.tmp')
            temporary.write_bytes(response.content)
            temporary.replace(path)
        table = pd.read_csv(path, sep='\t', usecols=['gene_id', 'TPM'], dtype={'gene_id':str}).iloc[1:]
        # Match the first-row removal in the original ENCODE downloader.
        if table.gene_id.duplicated().any():
            raise ValueError(f'Duplicate genes in {accession}')
        values = pd.to_numeric(table.TPM, errors='raise')
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f'Invalid TPM in {accession}')
        series.append(pd.Series(values.to_numpy(), index=table.gene_id, name=row.id))
        sources.append({'sample_id':row.id, 'accession':accession, 'url':url,
                        'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    raw = pd.concat(series, axis=1)
    annotation = pd.read_csv(annotations, dtype={'gene_id':str})
    if annotation.gene_id.duplicated().any():
        raise ValueError('Gene annotation must contain unique gene_id values')
    intervals = pd.concat([pd.read_csv(f'{intervals_prefix}{split}.{strand}.csv', sep='\t', dtype={'gene_id':str})
                          .assign(split=split)
                          for split in ['train','valid','test']
                          for strand in ['forward','reverse']], ignore_index=True)
    if intervals.gene_id.duplicated().any():
        raise ValueError('Interval splits must contain unique gene_id values')
    # Preserve the original coordinate-based annotation/interval intersection.
    genes = annotation.merge(intervals[['gene_id','chrom','TSS','TES','split']],
                             left_on=['gene_id','chromosome','TSS','TES'],
                             right_on=['gene_id','chrom','TSS','TES'], validate='1:1')
    genes = genes[genes.gene_id.isin(raw.index)].set_index('gene_id')
    raw = raw.reindex(genes.index)
    if raw.empty or raw.isna().any().any():
        raise ValueError('Empty matrix or missing TPM values; inspect sample gene sets')
    # All genes/samples participate before training-only low-expression filtering.
    normalized = qnorm.quantile_normalize(raw, ncpus=1)
    keep = (raw.max(axis=1) >= 2) | genes['split'].isin(['valid','test'])
    if (~keep).sum() > (genes['split'] != 'test').sum() / 5:
        keep[:] = True
    for name, frame in [('tpms',raw), ('qnorm_tpms',normalized)]:
        frame = frame.loc[keep].copy()
        frame['split'] = genes.loc[keep,'split']
        frame.to_csv(out / f'{organism}_{name}.csv', index_label='gene_id')
    mappings.to_csv(out / 'samples.csv', index=False)
    inputs = [Path(p) for p in file_mappings] + [Path(annotations)] + [
        Path(f'{intervals_prefix}{split}.{strand}.csv')
        for split in ['train','valid','test'] for strand in ['forward','reverse']]
    manifest = {'organism':organism, 'samples':len(mappings),
                'genes_normalized':len(raw), 'genes_saved':int(keep.sum()),
                'normalization':'qnorm.quantile_normalize on raw TPM, axis=1 (samples)',
                'source_commit':'6c10d3b7f2022488f40c6d3c52ddca338756cea5',
                'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                'encode_files':sources,
                'versions':{p:importlib.metadata.version(p) for p in ['qnorm','numpy','pandas']}}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Saved {len(mappings)} samples, {int(keep.sum())} genes to {out}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--file-mappings', nargs='+', required=True)
    parser.add_argument('--organism', choices=['human','mouse'], required=True)
    parser.add_argument('--annotations', required=True, help='Original genome gene annotation CSV')
    parser.add_argument('--intervals-prefix', required=True, help='Prefix before train.forward.csv, etc.')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    prepare(**vars(args))
