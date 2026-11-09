#!/usr/bin/env python3
"""Download ENCODE TPM targets; adapted from minjaf/ExpressionPredict.

Source: 6c10d3b7f2022488f40c6d3c52ddca338756cea5,
src/download_borzoi.py, src/qnorm_tpm_data.py and src/utils.py.
Uses recorded quantification accessions, so no BAM downloads are needed.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests


def download(file_mappings, organism, out):
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
        table = pd.read_csv(path, sep='\t', usecols=['gene_id', 'TPM'], dtype={'gene_id':str})
        prefix = 'ENSG' if organism == 'human' else 'ENSMUSG'
        table = table[table.gene_id.str.startswith(prefix)]
        if table.gene_id.duplicated().any():
            raise ValueError(f'Duplicate genes in {accession}')
        values = pd.to_numeric(table.TPM, errors='raise')
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f'Invalid TPM in {accession}')
        series.append(pd.Series(values.to_numpy(), index=table.gene_id, name=row.id))
        sources.append({'sample_id':row.id, 'accession':accession, 'url':url,
                        'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    raw = pd.concat(series, axis=1, join='inner')
    if raw.empty:
        raise ValueError('No common Ensembl genes across the selected samples')
    raw.to_csv(out / f'{organism}_tpms.csv', index_label='gene_id')
    mappings.to_csv(out / 'samples.csv', index=False)
    inputs = [Path(p) for p in file_mappings]
    manifest = {'organism':organism, 'samples':len(mappings), 'genes':len(raw),
                'gene_selection':'Common Ensembl gene IDs across all samples; no zero filling',
                'source_commit':'6c10d3b7f2022488f40c6d3c52ddca338756cea5',
                'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                'encode_files':sources}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Saved {len(mappings)} samples, {len(raw)} genes to {out}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--file-mappings', nargs='+', required=True)
    parser.add_argument('--organism', choices=['human','mouse'], required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    download(**vars(args))
