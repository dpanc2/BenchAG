#!/usr/bin/env python3
"""Compare 1 Mb AlphaGenome predictions with released RNA_SEQ TFRecord targets.

Requires TensorFlow for decoding, plus the standard benchmark dependencies.
The official metadata CSV includes padded channels in their original order.
No TPM conversion or experimental-sample regrouping is performed.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import struct
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from aggregate_gene_exons import read_exons, assign_windows
from evaluate_gene_expression import evaluate, plot_results


def now():
    return datetime.now(timezone.utc).isoformat()


def save_json(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def load_key(path):
    if os.environ.get('ALPHAGENOME_API_KEY'):
        return os.environ['ALPHAGENOME_API_KEY']
    for line in Path(path).read_text().splitlines():
        name, sep, value = line.strip().removeprefix('export ').partition('=')
        if sep and name.strip() == 'ALPHAGENOME_API_KEY':
            return value.strip().strip('\"\'')
    raise ValueError('API credential assignment not found')


def decode_record(raw):
    import tensorflow as tf
    example = tf.train.Example.FromString(raw)
    features = example.features.feature
    decoded = {}
    for key, dtype in [('dna_sequence', tf.float32), ('rna_seq', tf.bfloat16),
                       ('rna_seq_mask', tf.bool), ('rna_seq_strand', tf.int32)]:
        decoded[key] = tf.io.parse_tensor(features[key].bytes_list.value[0], dtype).numpy()
    decoded['chromosome'] = features['interval/chromosome'].bytes_list.value[0].decode()
    for key in ['start', 'end']:
        decoded[key] = features['interval/' + key].int64_list.value[0]
    return decoded


def read_exact(handle, size):
    chunks = []
    while size:
        chunk = handle.read(min(size, 4 * 1024 * 1024))
        if not chunk:
            raise EOFError('Truncated TFRecord')
        chunks.append(chunk)
        size -= len(chunk)
    return b''.join(chunks)


def stream_records(url):
    """Read a gzip TFRecord shard over HTTPS, checking native TFRecord CRCs."""
    import google_crc32c
    def checksum(data):
        crc = google_crc32c.value(data)
        return (((crc >> 15) | (crc << 17)) + 0xa282ead8) & 0xffffffff
    with urllib.request.urlopen(url, timeout=180) as response, gzip.GzipFile(fileobj=response) as handle:
        while True:
            size_bytes = handle.read(8)
            if not size_bytes:
                break
            if len(size_bytes) != 8:
                raise EOFError('Truncated TFRecord header')
            size_crc = struct.unpack('<I', read_exact(handle, 4))[0]
            if checksum(size_bytes) != size_crc:
                raise ValueError('TFRecord length checksum mismatch')
            raw = read_exact(handle, struct.unpack('<Q', size_bytes)[0])
            data_crc = struct.unpack('<I', read_exact(handle, 4))[0]
            if checksum(raw) != data_crc:
                raise ValueError('TFRecord data checksum mismatch')
            yield raw


def list_shards(prefix):
    entries, token = [], None
    while True:
        params = {'prefix': prefix, 'maxResults': 1000}
        if token:
            params['pageToken'] = token
        url = 'https://storage.googleapis.com/storage/v1/b/alphagenome-datasets/o?' + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=60) as response:
            page = json.load(response)
        entries.extend(page.get('items', []))
        token = page.get('nextPageToken')
        if not token:
            break
    entries = sorted([e for e in entries if e['name'].endswith('.gz.tfrecord')], key=lambda e: e['name'])
    if not entries:
        raise ValueError('No official target shards found')
    return entries


def validate_record(record, metadata):
    dna, targets = record['dna_sequence'], record['rna_seq']
    if dna.shape != (1048576, 4) or record['end'] - record['start'] != len(dna):
        raise ValueError('Expected an exact 1 Mb input window')
    if targets.ndim != 2 or targets.shape[1] != len(metadata):
        raise ValueError('Target channel count does not match official metadata')
    length = targets.shape[0]
    if length > len(dna) or (len(dna) - length) % 2:
        raise ValueError('Cannot center-align target and input windows')
    mask = np.broadcast_to(record['rna_seq_mask'], targets.shape)
    strands = np.asarray(record['rna_seq_strand']).reshape(-1)
    # Strand IDs in the released TFRecords: 0=unavailable, 1=+, 2=-, 3=unstranded.
    expected = metadata.strand.map({'+': 1, '-': 2, '.': 3}).to_numpy()
    available = mask.any(axis=0)
    if not np.array_equal(strands[available], expected[available]):
        raise ValueError('Official target strand IDs disagree with channel metadata')
    offset = (len(dna) - length) // 2
    return record['start'] + offset, record['start'] + offset + length, mask


def select_channels(metadata, strand):
    """Choose the gene-matched channel, or the single unstranded channel."""
    selected = metadata[(metadata.strand.isin([strand, '.'])) & (metadata.name.str.lower() != 'padding')].copy()
    # Some ontology/assay names have BOTH stranded and unstranded measurements.
    # They are separate target tracks; do not average them or overwrite either.
    duplicated = selected.name.duplicated(keep=False)
    selected['expression_name'] = selected['name']
    selected.loc[duplicated, 'expression_name'] = selected.loc[duplicated].apply(
        lambda row: row['name'] + (' | unstranded' if row['strand'] == '.' else ' | stranded'), axis=1)
    if selected.expression_name.duplicated().any():
        raise ValueError('Ambiguous channels even after separating strand protocols')
    return selected


def aggregate_record(record, metadata, assignments, prediction=None):
    start, end, mask = validate_record(record, metadata)
    target = record['rna_seq']
    observed, predicted = {}, {}
    if prediction is not None:
        if prediction.resolution != 1:
            raise ValueError('Need 1 bp RNA_SEQ predictions')
        interval = prediction.interval
        if (interval.chromosome, interval.start, interval.end) != (record['chromosome'], record['start'], record['end']):
            raise ValueError('Predicted and official input intervals differ')
        pred_meta = prediction.metadata.reset_index(drop=True)
        lookup = {}
        for i, row in pred_meta.iterrows():
            key = (row['name'], row['strand'])
            if key in lookup:
                raise ValueError('Duplicate prediction channel name/strand')
            lookup[key] = i
    for gid, assignment in assignments.items():
        obs_row, pred_row = {}, {}
        channels = select_channels(metadata, assignment['strand'])
        for index, channel in channels.iterrows():
            name = channel['expression_name']
            obs_sum, pred_sum, count = 0., 0., 0
            valid = True
            for a, b in assignment['regions']:
                if not (start <= a < b <= end):
                    raise ValueError('Exon outside target region')
                if not mask[a-start:b-start, index].all():
                    valid = False
                    break
                observed_values = target[a-start:b-start, index].astype(np.float64)
                if not np.isfinite(observed_values).all() or (observed_values < 0).any():
                    raise ValueError('Invalid observed coverage')
                obs_sum += observed_values.sum()
                count += b-a
                if prediction is not None:
                    key = (channel['name'], channel['strand'])
                    if key not in lookup:
                        raise ValueError(f'Missing predicted channel: {key}')
                    values = np.asarray(prediction.values[a-record['start']:b-record['start'], lookup[key]], dtype=float)
                    if not np.isfinite(values).all() or (values < 0).any():
                        raise ValueError('Invalid predicted coverage')
                    pred_sum += values.sum()
            obs_row[name] = obs_sum / count if valid and count else np.nan
            if prediction is not None:
                pred_row[name] = pred_sum / count if valid and count else np.nan
        observed[gid] = obs_row
        if prediction is not None:
            predicted[gid] = pred_row
    return observed, predicted


def sequence_from_onehot(dna):
    if not np.isfinite(dna).all() or not np.isin(dna, [0, 1]).all():
        raise ValueError('DNA must be finite one-hot or all-zero N rows')
    sums = dna.sum(axis=1)
    if not np.isin(sums, [0, 1]).all():
        raise ValueError('Invalid one-hot DNA')
    sequence = np.frombuffer(b'ACGT', dtype='S1')[dna.argmax(axis=1)].copy()
    sequence[sums == 0] = b'N'
    return sequence.tobytes().decode()


def predict_with_retries(model, sequence, **kwargs):
    """Retry only transient gRPC failures; checkpoints cover completed calls."""
    transient = {'UNAVAILABLE', 'RESOURCE_EXHAUSTED', 'DEADLINE_EXCEEDED'}
    for attempt in range(3):
        try:
            return model.predict_sequence(sequence, **kwargs)
        except Exception as error:
            code = error.code() if callable(getattr(error, 'code', None)) else None
            if getattr(code, 'name', None) not in transient or attempt == 2:
                raise
            print(f'Transient API status {code.name}; retry {attempt+1}/2', flush=True)
            time.sleep(5 * 2**attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument('--records', nargs='+', help='Local uncompressed or gzip TFRecords')
    sources.add_argument('--gcs-prefix', help='Exact official RNA_SEQ prefix')
    parser.add_argument('--metadata', required=True, help='Official padded RNA_SEQ metadata CSV, original order')
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--organism', choices=['HOMO_SAPIENS', 'MUS_MUSCULUS'], required=True)
    parser.add_argument('--model-version', default='FOLD_1')
    parser.add_argument('--api-env', required=True, help='Existing credential file; contents are never logged')
    parser.add_argument('--max-records', type=int, default=1, help='Bounded pilot by default; 0 means all records')
    parser.add_argument('--ontology-batch-size', type=int, default=16, help='Limit transferred output/memory per call')
    parser.add_argument('--export-every', type=int, default=25, help='Export cumulative matrices after this many gene-bearing windows')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    if args.max_records < 0 or args.ontology_batch_size < 1 or args.export_every < 1:
        parser.error('Invalid record/batch limit')
    if args.gcs_prefix and not args.gcs_prefix.startswith(f'v1/train/{args.model_version}/{args.organism}/TEST/RNA_SEQ/'):
        parser.error('Use the model-matched official TEST/RNA_SEQ prefix')
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    checkpoint = out / 'windows'
    checkpoint.mkdir(exist_ok=True)
    metadata = pd.read_csv(args.metadata, keep_default_na=False).reset_index(drop=True)
    if metadata.empty or not {'name', 'strand', 'ontology_curie'} <= set(metadata):
        raise ValueError('Invalid official metadata')
    identity = {'parameters': {k: v for k, v in vars(args).items() if k not in ['api_env', 'out']},
                'input_sha256': {str(Path(p).resolve()): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                                 for p in [args.metadata, args.gtf, __file__, str(Path(__file__).with_name('aggregate_gene_exons.py')),
                                           str(Path(__file__).with_name('evaluate_gene_expression.py'))] + (args.records or [])}}
    manifest_path = out / 'run.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text())['identity'] != identity:
        raise ValueError('Run inputs or code changed; use a new output directory')
    manifest = {'identity': identity, 'started_utc': now(), 'status': 'running', 'pid': os.getpid(),
                'checkpoint_windows': 0, 'api_calls_this_process': 0, 'paper_reproduction_certified': False,
                'channel_aggregation': 'logical track: matching gene strand or unstranded channel',
                'versions': {p: importlib.metadata.version(p) for p in ['alphagenome', 'numpy', 'tensorflow-cpu']}}
    save_json(manifest_path, manifest)
    # Prevent proxies from interfering with gRPC; HTTPS data download remains direct.
    for name in ['HTTP_PROXY', 'HTTPS_PROXY', 'ALL_PROXY', 'http_proxy', 'https_proxy', 'all_proxy']:
        os.environ.pop(name, None)
    import tensorflow as tf
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    model = dna_client.create(load_key(args.api_env), model_version=getattr(dna_client.ModelVersion, args.model_version), timeout=600)
    organism = getattr(dna_client.Organism, args.organism)
    api_metadata = model.output_metadata(organism=organism).rna_seq
    target_pairs = set(zip(metadata.loc[metadata.name.str.lower() != 'padding', 'name'], metadata.loc[metadata.name.str.lower() != 'padding', 'strand']))
    if target_pairs != set(zip(api_metadata.name, api_metadata.strand)):
        raise ValueError('API channel set differs from the official target metadata')
    api_metadata.to_csv(out / 'api_metadata.csv', index=False)
    genes = read_exons(args.gtf)
    genes_by_chrom = {}
    for gid, gene in genes.items():
        genes_by_chrom.setdefault(gene['chrom'], {})[gid] = gene
    if args.records:
        def iterator():
            for path in args.records:
                compression = 'GZIP' if path.endswith('.gz.tfrecord') or path.endswith('.gz') else ''
                for i, raw in enumerate(tf.data.TFRecordDataset([path], compression_type=compression)):
                    yield path, i, raw.numpy()
    else:
        inventory_path = out / 'target_objects.json'
        if inventory_path.exists():
            objects = json.loads(inventory_path.read_text())
        else:
            objects = list_shards(args.gcs_prefix)
            save_json(inventory_path, objects)
        def iterator():
            for obj in objects:
                url = 'https://storage.googleapis.com/alphagenome-datasets/' + obj['name'] + '?generation=' + obj['generation']
                for i, raw in enumerate(stream_records(url)):
                    yield obj['name'], i, raw
    def export():
        # Resolve duplicated genes across overlapping windows deterministically.
        chosen = {}
        for path in sorted(checkpoint.glob('*.json')):
            if '.batch' in path.name:
                continue
            window = json.loads(path.read_text())
            for gid, gene in window['genes'].items():
                if gid not in chosen or gene['exon_bases'] > chosen[gid][0]['exon_bases']:
                    chosen[gid] = (gene, window['observed'][gid], window['predicted'][gid])
        if not chosen:
            manifest['evaluation_status'] = 'no_gene_predictions'
            return
        pred = pd.DataFrame.from_dict({g: v[2] for g, v in chosen.items()}, orient='index').rename_axis('gene_id')
        obs = pd.DataFrame.from_dict({g: v[1] for g, v in chosen.items()}, orient='index').rename_axis('gene_id')
        pred.to_csv(out / 'predictions.csv')
        obs.to_csv(out / 'observed.csv')
        manifest['genes_saved'] = len(chosen)
        try:
            summary, scores, coverage = evaluate(np.log1p(pred), np.log1p(obs))
        except ValueError as error:
            manifest['evaluation_status'] = 'not_ready'
            manifest['evaluation_pending_reason'] = str(error)
            save_json(manifest_path, manifest)
            return
        summary.to_csv(out / 'summary.csv', index=False)
        scores.to_csv(out / 'correlations.csv', index=False)
        save_json(out / 'coverage.json', coverage)
        plot_results(scores, out / 'gene_expression.png', f'{args.organism}: official targets, 1 Mb, {args.model_version}')
        manifest['evaluation_status'] = 'available'
        manifest.pop('evaluation_pending_reason', None)
        print(summary.to_string(index=False), flush=True)
    try:
        processed = 0
        for source, record_index, raw in iterator():
            manifest['records_scanned'] = manifest.get('records_scanned', 0) + 1
            tag = hashlib.sha256(f'{source}:{record_index}'.encode()).hexdigest()[:20]
            window_path = checkpoint / (tag + '.json')
            if not window_path.exists():
                record = decode_record(raw)
                start, end, mask = validate_record(record, metadata)
                region = (record['chromosome'], start, end)
                assignments = assign_windows(genes_by_chrom.get(record['chromosome'], {}), [region])
                print(f'Window {source}:{record_index} {region}; {len(assignments)} eligible genes', flush=True)
                if not assignments:
                    save_json(window_path, {'source': source, 'record_index': record_index,
                                            'input_interval': [record['chromosome'], record['start'], record['end']],
                                            'target_interval': list(region), 'genes': {}, 'observed': {}, 'predicted': {}})
                    manifest['skipped_gene_empty_windows'] = manifest.get('skipped_gene_empty_windows', 0) + 1
                    save_json(manifest_path, manifest)
                    continue
                observed, _ = aggregate_record(record, metadata, assignments)
                predicted = {g: {} for g in assignments}
                available = mask.any(axis=0) & (metadata.name.str.lower() != 'padding').to_numpy()
                ontologies = list(dict.fromkeys(metadata.loc[available, 'ontology_curie']))
                if not all(ontologies):
                    raise ValueError('Available target channel without ontology')
                sequence = sequence_from_onehot(record.pop('dna_sequence'))
                # Keep a lightweight shape proxy after extracting DNA for validation.
                record['dna_sequence'] = np.broadcast_to(np.uint8(0), (1048576, 4))
                for batch_start in range(0, len(ontologies), args.ontology_batch_size):
                    batch = ontologies[batch_start:batch_start + args.ontology_batch_size]
                    cache_path = checkpoint / f'{tag}.batch{batch_start:04d}.json'
                    if cache_path.exists():
                        partial = json.loads(cache_path.read_text())
                    else:
                        output = predict_with_retries(model, sequence, interval=genome.Interval(chromosome=record['chromosome'], start=record['start'], end=record['end']),
                                                        organism=organism, requested_outputs=[dna_client.OutputType.RNA_SEQ], ontology_terms=batch)
                        manifest['api_calls_this_process'] += 1
                        subset = metadata.copy()
                        subset.loc[~subset.ontology_curie.isin(batch), 'name'] = 'Padding'
                        _, partial = aggregate_record(record, subset, assignments, output.rna_seq)
                        save_json(cache_path, partial)
                        del output
                    for gid, values in partial.items():
                        predicted[gid].update(values)
                    manifest['current_window'] = {'source': source, 'record': record_index, 'completed_ontology_batches': batch_start // args.ontology_batch_size + 1,
                                                  'total_ontology_batches': (len(ontologies) + args.ontology_batch_size - 1) // args.ontology_batch_size}
                    save_json(manifest_path, manifest)
                    print(f'  completed batch {batch_start // args.ontology_batch_size + 1}', flush=True)
                save_json(window_path, {'source': source, 'record_index': record_index, 'input_interval': [record['chromosome'], record['start'], record['end']],
                                        'target_interval': list(region), 'target_record_sha256': hashlib.sha256(raw).hexdigest(),
                                        'genes': assignments, 'observed': observed, 'predicted': predicted})
            if not json.loads(window_path.read_text())['genes']:
                continue
            processed += 1
            manifest['gene_bearing_windows'] = processed
            manifest['checkpoint_windows'] = len(list(checkpoint.glob('*.json'))) - len(list(checkpoint.glob('*.batch*.json')))
            save_json(manifest_path, manifest)
            if processed % args.export_every == 0:
                export()
            if args.max_records and processed >= args.max_records:
                break
        export()
        manifest.update(status='completed' if manifest.get('evaluation_status') == 'available' else 'completed_without_metrics', finished_utc=now())
        save_json(manifest_path, manifest)
    except Exception as error:
        manifest.update(status='failed', error_type=type(error).__name__, error=str(error), finished_utc=now())
        save_json(manifest_path, manifest)
        raise


if __name__ == '__main__':
    main()
