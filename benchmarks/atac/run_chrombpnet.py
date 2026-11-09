#!/usr/bin/env python3
"""Run the colleague's ChromBPNet prediction and metrics with saved region profiles."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyBigWig
import pysam

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/atac'
BED = DATA / 'source/alpha/fold1_test_regions.bed'
TRUTH = DATA / 'ENCFF650PHA.bigWig'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genome', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--limit', type=int)
    args = p.parse_args();out = args.out.resolve();out.mkdir(parents=True, exist_ok=True)
    model = DATA / 'models/chrombpnet_recompiled.h5'
    identity = {'source_commit':'9e6a91157702b80efbc3b43f23c06e21f7bae856',
                'model_sha256':hashlib.sha256(model.read_bytes()).hexdigest(),
                'bed_sha256':hashlib.sha256(BED.read_bytes()).hexdigest(),
                'genome':str(args.genome.resolve()),'limit':args.limit,
                'prediction_center':'summit','evaluation_center':'peak midpoint',
                'input_length':2114,'output_length':1000,'reverse_complement':True}
    manifest = out / 'inputs.json'
    if manifest.exists() and json.loads(manifest.read_text()) != identity:
        raise ValueError('Inputs/settings changed; use a new output directory')
    manifest.write_text(json.dumps(identity, indent=2)+'\n')
    if not Path(str(args.genome)+'.fai').exists():
        pysam.faidx(str(args.genome))
    genome = pysam.FastaFile(str(args.genome))
    sizes = out / 'chrom.sizes'
    sizes.write_text(''.join(f'{c}\t{n}\n' for c,n in zip(genome.references,genome.lengths)))
    genome.close()
    selected = BED
    if args.limit:
        selected = out / 'selected_regions.bed'
        pd.read_csv(BED, sep='\t',header=None).head(args.limit).to_csv(selected,sep='\t',header=False,index=False)
    predictor = load('colleague_chrom', DATA/'source/chrombpnet/chrom_bw.py')
    original = predictor.predict_averaged_profile
    cache = out/'regions';cache.mkdir(exist_ok=True)

    def cached(model, genome, chrom, start, end, inputlen, outputlen):
        path = cache/f'{chrom}_{start}_{end}.npy'
        if path.exists():
            return np.load(path)
        values = original(model, genome, chrom, start, end, inputlen, outputlen)
        if values.shape != (1000,) or not np.isfinite(values).all():
            raise ValueError('Invalid ChromBPNet prediction')
        temporary = path.with_suffix('.tmp.npy');np.save(temporary,values);temporary.replace(path)
        return values

    predictor.predict_averaged_profile = cached
    predicted_bw = out/'chrombpnet_averaged.bw'
    completion = out/'bigwig.complete'
    if not (predicted_bw.exists() and completion.exists()):
        sys.argv = ['chrom_bw.py','--model',str(model),'--genome',str(args.genome),
                    '--regions',str(selected),'--output_prefix',str(out/'chrombpnet'),
                    '--inputlen','2114','--outputlen','1000','--chrom_sizes',str(sizes)]
        predictor.main()
        missing = []
        for row in pd.read_csv(selected, sep='\t', header=None).itertuples(index=False, name=None):
            center = int(row[1]) + int(row[9])
            if not (cache/f'{row[0]}_{center-1057}_{center+1057}.npy').exists():
                missing.append([row[0],int(row[1]),int(row[2])])
        if missing:
            (out/'missing_predictions.json').write_text(json.dumps(missing,indent=2)+'\n')
            raise RuntimeError('Some ChromBPNet regions lack valid predictions; rerun to retry')
        if not predicted_bw.exists():
            raise RuntimeError('Predicted BigWig was not produced')
        completion.write_text('All requested region profiles saved and BigWig generation finished.\n')
    metric = load('colleague_chrom_metrics',DATA/'source/chrombpnet/chrom_metrics.py')
    metric.batching = False
    result = metric.evaluate_sample('ATAC_IMR-90',str(TRUTH),str(predicted_bw),str(selected),
                                    500,1,str(out),1000,1e-16)
    missing=0
    with pyBigWig.open(str(predicted_bw)) as bw:
        for row in pd.read_csv(selected,sep='\t',header=None).itertuples(index=False,name=None):
            center=(int(row[1])+int(row[2]))//2
            missing += int(np.isnan(bw.values(row[0],center-500,center+500,numpy=True)).sum())
    summary={k:v for k,v in result.items() if not k.startswith('per_sample_')}
    summary.update(model='ChromBPNet fold0',missing_prediction_base_slots=missing,
                   status='complete' if result['n_failed_regions']==0 else 'incomplete')
    (out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame([summary]).to_csv(out/'summary.csv',index=False)
    if result['n_failed_regions']:
        raise RuntimeError('Some evaluation regions failed')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
