#!/usr/bin/env python3
"""Aggregate matched observed/predicted RNA-seq bigWigs over GTF exon unions.

Optional --api obtains predicted coverage from AlphaGenome instead of bigWigs.
All coordinates internally are zero-based half-open; GTF is converted on read.
"""
import argparse
from contextlib import ExitStack
import gzip
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd


def merge_exons(exons):
    merged = []
    for start, end in sorted(set(exons)):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def read_exons(path):
    genes = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'exon':
                continue
            match = re.search(r'(?:^|;\s*)gene_id "([^"]+)"', fields[8])
            if not match:
                raise ValueError('Exon row lacks gene_id')
            gid = match[1]
            chrom, strand = fields[0], fields[6]
            start, end = int(fields[3]) - 1, int(fields[4])
            if strand not in ['+', '-'] or start < 0 or end <= start:
                raise ValueError('Invalid exon coordinate/strand')
            gene = genes.setdefault(gid, {'chrom': chrom, 'strand': strand, 'exons': set()})
            if (chrom, strand) != (gene['chrom'], gene['strand']):
                raise ValueError(f'{gid}: inconsistent chromosome/strand')
            gene['exons'].add((start, end))
    if not genes:
        raise ValueError('Annotation contains no exons')
    return genes


def assign_windows(genes, windows, min_fraction=.5):
    """Count distinct annotated exon intervals fully contained in each window.

    Select the qualifying window covering most union-exon bases, with BED order
    breaking ties. This explicit choice avoids counting a gene twice.
    """
    assigned = {}
    by_chrom = {}
    for i, (chrom, start, end) in enumerate(windows):
        by_chrom.setdefault(chrom, []).append((i, start, end))
    for gid, gene in genes.items():
        union = merge_exons(gene['exons'])
        best = None
        for i, start, end in by_chrom.get(gene['chrom'], []):
            contained = sum(start <= a and b <= end for a, b in gene['exons'])
            fraction = contained / len(gene['exons'])
            if fraction < min_fraction:
                continue
            clipped = [(max(a, start), min(b, end)) for a, b in union if a < end and b > start]
            bases = sum(b-a for a, b in clipped)
            if bases and (best is None or bases > best['exon_bases']):
                best = {'gene_id': gid, 'window_index': i, 'exon_fraction': fraction,
                        'exon_bases': bases, 'regions': clipped, 'chrom': gene['chrom'], 'strand': gene['strand']}
        if best is not None:
            assigned[gid] = best
    return assigned


def mean_coverage(reader, chrom, regions, missing='error'):
    """Base-weighted mean, with no double counting of overlapping isoform exons."""
    total, count = 0., 0
    for start, end in regions:
        values = np.asarray(reader.values(chrom, start, end, numpy=True), dtype=float)
        if len(values) != end-start:
            raise ValueError('Coverage length does not match exon')
        if np.isinf(values).any() or (values[np.isfinite(values)] < 0).any():
            raise ValueError('RNA coverage must be finite and nonnegative')
        if np.isnan(values).any():
            if missing == 'error':
                raise ValueError('Missing bigWig bases; use --missing-bases zero only for documented sparse tracks')
            values = np.nan_to_num(values, nan=0.)
        total += values.sum()
        count += len(values)
    return total / count


def api_gene_means(model, window, assignments, tracks, organism):
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    chrom, start, end = window
    output = model.predict_interval(
        genome.Interval(chromosome=chrom, start=start, end=end),
        organism=organism, requested_outputs=[dna_client.OutputType.RNA_SEQ], ontology_terms=None)
    data = output.rna_seq
    if data is None or data.resolution != 1:
        raise ValueError('Need 1 bp RNA_SEQ predictions')
    if (data.interval.chromosome, data.interval.start, data.interval.end) != window:
        raise ValueError('API returned a different output interval')
    result = {}
    for gid, gene in assignments.items():
        row = {}
        for track in tracks.itertuples(index=False):
            if track.strand != gene['strand']:
                continue
            mask = (data.metadata['name'] == track.track_id) & (data.metadata['strand'] == track.strand)
            indices = np.flatnonzero(mask)
            if len(indices) != 1:
                raise ValueError(f'Expected one API track: {track.track_id} {track.strand}; got {len(indices)}')
            total = 0.
            for a, b in gene['regions']:
                values = np.asarray(data.values[a-start:b-start, indices[0]], dtype=float)
                if not np.isfinite(values).all() or (values < 0).any():
                    raise ValueError('Invalid API coverage values')
                total += values.sum()
            row[track.track_id] = total / gene['exon_bases']
        result[gid] = row
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gtf', required=True)
    p.add_argument('--windows', required=True, help='BED3: evaluation intervals, not TSS positions')
    p.add_argument('--tracks', required=True, help='CSV: track_id,strand,observed_bigwig[,predicted_bigwig]')
    p.add_argument('--min-exon-fraction', type=float, default=.5)
    p.add_argument('--missing-bases', choices=['error', 'zero'], default='error')
    p.add_argument('--api', action='store_true', help='Use API RNA_SEQ output; consumes quota')
    p.add_argument('--model-version', default='FOLD_1')
    p.add_argument('--organism', choices=['human', 'mouse'], default='human')
    p.add_argument('--assembly', required=True)
    p.add_argument('--annotation-release', required=True)
    p.add_argument('--split-provenance', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    if not 0 < args.min_exon_fraction <= 1:
        p.error('min-exon-fraction must be in (0,1]')
    windows = pd.read_csv(args.windows, sep='\t', comment='#', header=None).iloc[:, :3]
    windows.columns = ['chrom', 'start', 'end']
    if windows.empty or windows.isna().any().any():
        p.error('Empty/invalid BED')
    for col in ['start', 'end']:
        if not np.equal(windows[col], windows[col].astype(int)).all():
            p.error('BED coordinates must be integers')
    if (windows.start < 0).any() or (windows.end <= windows.start).any():
        p.error('Invalid BED coordinates')
    windows = list(windows.itertuples(index=False, name=None))
    tracks = pd.read_csv(args.tracks)
    required = {'track_id', 'strand', 'observed_bigwig'} | (set() if args.api else {'predicted_bigwig'})
    if not required <= set(tracks) or tracks[list(required)].isna().any().any() or tracks.empty:
        p.error('Invalid tracks manifest')
    if tracks.duplicated(['track_id', 'strand']).any() or not tracks.strand.isin(['+', '-']).all():
        p.error('Track/strand pairs must be unique and stranded')
    # A logical track is represented by its matching-strand channel for each gene.
    if any(set(group.strand) != {'+', '-'} for _, group in tracks.groupby('track_id')):
        p.error('Each logical track needs both + and - rows')
    assigned = assign_windows(read_exons(args.gtf), windows, args.min_exon_fraction)
    if not assigned:
        p.error('No genes satisfy the exon coverage threshold')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        p.error('Use an empty output directory')
    model = organism = None
    if args.api:
        if not os.environ.get('ALPHAGENOME_API_KEY'):
            p.error('Set ALPHAGENOME_API_KEY')
        from alphagenome.models import dna_client
        if any(end-start not in [16384, 131072, 524288, 1048576] for _, start, end in windows):
            p.error('API windows must have a supported sequence length')
        organism = getattr(dna_client.Organism, 'HOMO_SAPIENS' if args.organism == 'human' else 'MUS_MUSCULUS')
        model = dna_client.create(os.environ['ALPHAGENOME_API_KEY'],
                                  model_version=getattr(dna_client.ModelVersion, args.model_version))
    import pyBigWig
    observed, predicted = {}, {}
    # Open each file once; repeated +/- manifest rows can share a file.
    with ExitStack() as stack:
        handles = {}
        def reader(path):
            if path not in handles:
                handle = pyBigWig.open(path)
                if handle is None:
                    raise ValueError(f'Cannot open bigWig: {path}')
                handles[path] = stack.enter_context(handle)
            return handles[path]
        for gid, gene in assigned.items():
            observed[gid], predicted[gid] = {}, {}
            for track in tracks.itertuples(index=False):
                if track.strand == gene['strand']:
                    observed[gid][track.track_id] = mean_coverage(reader(track.observed_bigwig), gene['chrom'], gene['regions'], args.missing_bases)
                    if not args.api:
                        predicted[gid][track.track_id] = mean_coverage(reader(track.predicted_bigwig), gene['chrom'], gene['regions'], args.missing_bases)
        if args.api:
            for i, window in enumerate(windows):
                subset = {g: a for g, a in assigned.items() if a['window_index'] == i}
                if subset:
                    predicted.update(api_gene_means(model, window, subset, tracks, organism))
                    print(f'Predicted window {i+1}/{len(windows)}: {len(subset)} genes', flush=True)
                    # Save completed windows even if a later API call fails.
                    pd.DataFrame.from_dict(predicted, orient='index').rename_axis('gene_id').to_csv(out / 'predictions.partial.csv')
    for name, values in [('predictions', predicted), ('observed', observed)]:
        pd.DataFrame.from_dict(values, orient='index').rename_axis('gene_id').to_csv(out / f'{name}.csv')
    pd.DataFrame([{k: v for k, v in a.items() if k != 'regions'} for a in assigned.values()]).to_csv(out / 'genes.csv', index=False)
    files = [args.gtf, args.windows, args.tracks, __file__]
    coverage_paths = sorted(set(tracks.observed_bigwig) | (set() if args.api else set(tracks.predicted_bigwig)))
    # Hash annotation/configs; bigWigs can be enormous or remote, record paths/sizes.
    manifest = {'arguments': vars(args), 'n_genes': len(assigned), 'n_tracks': tracks.track_id.nunique(),
                'exon_definition': 'union of all annotated exon bases; unique exon intervals counted for 50% filter',
                'deduplication': 'maximum covered exon bases, then BED order',
                'output_units': 'linear mean read coverage; apply log1p in evaluator',
                'paper_reproduction_certified': False,
                'sha256': {str(Path(f).resolve()): hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in files},
                'coverage_files': [{'path': f, 'size': Path(f).stat().st_size if Path(f).is_file() else None} for f in coverage_paths]}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Saved {len(assigned)} genes to {out}')


if __name__ == '__main__':
    main()
