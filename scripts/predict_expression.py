#!/usr/bin/env python3
"""Predict expression from a TSS table and sample mapping; ground truth is optional."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
WINDOWS = {'16kb': (16384, 10001), '128kb': (131072, 100001), '1mb': (1048576, 200001)}

def read_genes(path, organism):
    import pandas as pd
    table = pd.read_csv(path, sep=None, engine='python', dtype={'gene_id': str})
    required = {'gene_id', 'chrom', 'TSS', 'gene_strand'}
    if not required <= set(table):
        raise ValueError(f'Gene table needs columns: {", ".join(sorted(required))}')
    if table.empty or table[list(required)].isna().any().any():
        raise ValueError('Gene table is empty or has missing required values')
    table['gene_id_unversioned'] = table.gene_id.str.split('.').str[0]
    if table.gene_id_unversioned.duplicated().any():
        raise ValueError('Use one TSS per gene; duplicate unversioned gene IDs found')
    if not table.gene_strand.isin(['+', '-']).all():
        raise ValueError('gene_strand must be + or -')
    prefix = 'ENSG' if organism == 'human' else 'ENSMUSG'
    if not table.gene_id_unversioned.str.startswith(prefix).all():
        raise ValueError('Gene IDs do not match the selected organism')
    table['TSS'] = pd.to_numeric(table.TSS, errors='raise')
    if (table.TSS < 0).any() or (table.TSS % 1 != 0).any():
        raise ValueError('TSS must contain nonnegative integer, zero-based coordinates')
    return table

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genes', required=True, help='CSV/TSV: gene_id,chrom,TSS,gene_strand')
    parser.add_argument('--mapping', required=True, help='track_to_first_id.json from prepare_expression_mapping.py')
    parser.add_argument('--organism', choices=['human', 'mouse'], default='human')
    parser.add_argument('--window', choices=list(WINDOWS), default='128kb')
    parser.add_argument('--model-version', default='FOLD_1')
    parser.add_argument('--out', required=True, help='Output directory; rerun unchanged command to resume')
    parser.add_argument('--check', action='store_true', help='Validate without API calls')
    args = parser.parse_args()
    table = read_genes(args.genes, args.organism)
    mapping = json.loads(Path(args.mapping).read_text())
    if not mapping or not all(isinstance(k, str) and isinstance(v, str) for k, v in mapping.items()):
        parser.error('Mapping must contain ontology-assay names mapped to representative sample IDs')
    if len(set(mapping.values())) != len(mapping):
        parser.error('Representative sample IDs must be unique across output names')
    if args.check:
        print(f'Inputs valid: {len(table)} genes, {len(mapping)} ontology-assay names. No API calls made.')
        return
    if not os.environ.get('ALPHAGENOME_API_KEY'):
        parser.error('Set ALPHAGENOME_API_KEY before running predictions')
    folder = ROOT / 'benchmarks/expression/scripts'
    source = folder / 'Alphagenome_prediction_benchmarking_json.py'
    identity = {'organism': args.organism, 'window': args.window, 'model_version': args.model_version,
                'client_version': importlib.metadata.version('alphagenome'),
                'sha256': {str(Path(p).resolve()): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [args.genes, args.mapping, source, __file__]}}
    out = Path(args.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / 'prediction_manifest.json'
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if previous['identity'] != identity:
            parser.error('Inputs or settings changed. Use another --out directory to keep checkpoints separate.')
    manifest = {'identity': identity, 'started_utc': datetime.now(timezone.utc).isoformat(), 'status': 'running'}
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    os.environ['ALPHAGENOME_ONTOLOGY_JSON'] = str(Path(args.mapping).resolve())
    sys.path.insert(0, str(folder))
    spec = importlib.util.spec_from_file_location('expression_runner', source)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.SEQUENCE_LENGTH, module.GENE_MASK_WIDTH = WINDOWS[args.window]
    module.ORGANISM = getattr(module.dna_client.Organism, 'HOMO_SAPIENS' if args.organism == 'human' else 'MUS_MUSCULUS')
    intervals = module.build_gene_intervals_dict(table[table.gene_strand == '+'], table[table.gene_strand == '-'])
    try:
        matrix = module.get_fold_intervals_for_intervals(args.model_version, os.environ['ALPHAGENOME_API_KEY'], intervals, None,
                            checkpoint_path=str(out / 'checkpoint.csv'), checkpoint_every=5)
    except Exception as error:
        manifest.update(status='failed', error_type=type(error).__name__, finished_utc=datetime.now(timezone.utc).isoformat())
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        raise
    # Public outputs use AlphaGenome names; internal checkpoint layout remains unchanged.
    matrix = matrix.rename(index={value: name for name, value in mapping.items()})
    frame = matrix.T
    versioned = dict(zip(table.gene_id_unversioned, table.gene_id))
    frame.index = frame.index.map(versioned); frame.index.name = 'gene_id'
    frame.to_csv(out / 'predictions.csv')
    manifest.update(status='completed', genes_requested=len(table), genes_saved=len(frame),
                    ontology_names_requested=len(mapping), finished_utc=datetime.now(timezone.utc).isoformat())
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Saved {len(frame)} genes to {out / "predictions.csv"}')

if __name__ == '__main__':
    try: main()
    except (ValueError, FileNotFoundError) as exc: raise SystemExit(f'Input error: {exc}')
