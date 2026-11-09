#!/usr/bin/env python3
"""Run the colleague's ATAC evaluation with explicit folds and resumable predictions."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import types
import time
import threading

import numpy as np
import pandas as pd
import pyBigWig
from alphagenome.models import dna_client

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'data/atac/source/alpha/alpha.py'
BED = ROOT / 'data/atac/source/alpha/fold1_test_regions.bed'
TRUTH = ROOT / 'data/atac/ENCFF650PHA.bigWig'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-version', choices=['FOLD_0', 'FOLD_1'], required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--limit', type=int, help='First N regions, for a smoke test only')
    parser.add_argument('--workers', type=int, choices=range(1, 9), default=2)
    args = parser.parse_args()
    out = args.out.resolve();out.mkdir(parents=True, exist_ok=True)
    cache = out / 'regions';cache.mkdir(exist_ok=True)
    identity = {'model_version':args.model_version, 'input_length':1048576,
                'eval_window':1000, 'ontology':'EFO:0001196', 'output':'ATAC',
                'reverse_complement':'upstream forward/reverse averaging',
                'limit':args.limit, 'source_commit':'9e6a91157702b80efbc3b43f23c06e21f7bae856',
                'input_sha256':{str(p):digest(p) for p in [SOURCE, BED, TRUTH]}}
    manifest = out / 'inputs.json'
    if manifest.exists() and json.loads(manifest.read_text()) != identity:
        raise ValueError('Inputs/settings changed; use a new output directory')
    manifest.write_text(json.dumps(identity, indent=2)+'\n')
    # torch is imported but unused in the original API script.
    upstream = types.ModuleType('colleague_atac')
    code = SOURCE.read_text().replace('import torch\n', '')
    exec(compile(code, str(SOURCE), 'exec'), upstream.__dict__)
    original_predict = upstream.predict_with_reverse_complement

    def cached_predict(model, pred_interval, eval_interval, biosample, data_type):
        name = f'{eval_interval.chromosome}_{eval_interval.start}_{eval_interval.end}.npz'
        path = cache / name
        if path.exists():
            with np.load(path) as saved:
                return saved['prediction']
        for attempt in range(3):
            try:
                values = np.asarray(original_predict(model, pred_interval, eval_interval, biosample, data_type))
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 ** (attempt + 1))
        if values.shape != (1000,) or not np.isfinite(values).all():
            raise ValueError(f'Unexpected prediction shape/values: {values.shape}')
        with pyBigWig.open(str(TRUTH)) as truth_bw:
            observed = np.nan_to_num(truth_bw.values(eval_interval.chromosome, eval_interval.start,
                                                    eval_interval.end, numpy=True))
        temporary = path.with_suffix(f'.{threading.get_ident()}.tmp.npz')
        np.savez_compressed(temporary, prediction=values, truth=observed)
        temporary.replace(path)
        return values

    upstream.predict_with_reverse_complement = cached_predict
    upstream.batching = False
    selected = BED
    table = pd.read_csv(BED, sep='\t', header=None)
    if args.limit:
        table = table.head(args.limit)
        selected = out / 'selected_regions.bed'
        table.to_csv(selected, sep='\t', header=False, index=False)
    table.to_csv(out/'evaluation_regions.bed', sep='\t', header=False, index=False)
    key = os.environ.get('ALPHAGENOME_API_KEY')
    if not key:
        raise ValueError('Set ALPHAGENOME_API_KEY before running')
    model = dna_client.create(key, model_version=getattr(dna_client.ModelVersion, args.model_version))
    # Prefetch unique windows; the original evaluation then reads the saved arrays.
    windows = {}
    for row in table.itertuples(index=False, name=None):
        interval = upstream.genome.Interval(str(row[0]), int(row[1]), int(row[2]))
        evaluation = upstream.get_evaluation_window(interval, 1000)
        windows[(evaluation.chromosome, evaluation.start, evaluation.end)] = evaluation
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        tasks = {pool.submit(cached_predict, model, w.resize(1048576), w,
                             'EFO:0001196', 'ATAC'):key for key,w in windows.items()}
        for n,future in enumerate(as_completed(tasks), 1):
            try:
                future.result()
            except Exception as error:
                failures.append({'window':tasks[future], 'error':str(error)})
            if n % 25 == 0 or n == len(tasks):
                print(f'Saved/checked {n}/{len(tasks)} windows; failures: {len(failures)}', flush=True)
    if failures:
        (out/'errors.json').write_text(json.dumps(failures, indent=2)+'\n')
        raise RuntimeError('Prediction failures recorded; rerun the same command to retry missing windows')
    result = upstream.evaluate_sample('ATAC_IMR-90', 'EFO:0001196', str(TRUTH), str(selected),
                                     'ATAC', model, 1048576, 500, 1, str(out), 1000, 1e-16)
    summary = {k:v for k,v in result.items() if not k.startswith('per_sample_')}
    summary.update(model_version=args.model_version,
                   status='complete' if result['n_failed_regions']==0 else 'incomplete',
                   held_out=(args.model_version=='FOLD_1'))
    (out/'metrics.json').write_text(json.dumps(summary, indent=2)+'\n')
    pd.DataFrame([summary]).to_csv(out/'summary.csv', index=False)
    if result['n_failed_regions']:
        raise RuntimeError('Some regions failed; saved predictions can be resumed with the same command')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
