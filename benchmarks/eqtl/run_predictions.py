#!/usr/bin/env python3
"""Run the colleague's prediction code with per-variant caching for safe resume."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time
import grpc
import pandas as pd
from alphagenome.models import variant_scorers


def run(variants, mapping, model, out, source, cache, limit):
    variants, mapping, out, source, cache = map(Path,[variants,mapping,out,source,cache])
    out.mkdir(parents=True,exist_ok=True); cache=cache/model; cache.mkdir(parents=True,exist_ok=True)
    code=(source/'code/make_preds.py').read_text()
    inputs=pd.read_csv(variants,sep='\t')
    if limit:
        inputs=inputs[inputs.variant.isin(inputs.variant.drop_duplicates().head(limit))]
    target=out/variants.name; inputs.to_csv(target,sep='\t',index=False)
    all_genes={}
    for name in ['snv_eqtl','indel_credible_sets_split','indel_credible_sets_merged']:
        f=pd.read_csv(variants.parent/f'{name}.tsv',sep='\t')
        for v,g in f.groupby('variant').gene_id:
            all_genes.setdefault(v,set()).update(g)
    manifest={'model':model,'input_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
              'mapping_sha256':hashlib.sha256(mapping.read_bytes()).hexdigest(),
              'source_sha256':hashlib.sha256(code.encode()).hexdigest(),
              'sdk':importlib.metadata.version('alphagenome'),'cache':'Fresh API scores only; no colleague predictions loaded',
              'changes':'Only API/tidy block wrapped for retries and caching; original aggregation and output code retained'}
    old=out/'inputs.json'
    if old.exists() and json.loads(old.read_text())!=manifest:
        raise ValueError('Input/settings changed; choose a new output folder')
    old.write_text(json.dumps(manifest,indent=2)+'\n')
    def cached_tidy(client,interval,variant,scorer,variant_id):
        genes=all_genes[variant_id]
        identity=json.dumps([variant_id,sorted(genes),manifest['sdk'],manifest['source_sha256']],sort_keys=True)
        path=cache/(hashlib.sha256(identity.encode()).hexdigest()+'.parquet')
        if path.exists(): return pd.read_parquet(path)
        for attempt in range(10):
            try:
                score=client.score_variant(interval,variant,[scorer])[0]
                break
            except grpc.RpcError:
                if attempt==9:raise
                time.sleep(60)
        tidy=variant_scorers.tidy_scores([[score]],match_gene_strand=True,include_extended_metadata=True)
        tidy=tidy.loc[tidy.gene_id.str.split('.').str[0].isin(genes), ['gene_id','gtex_tissue','raw_score']].copy()
        temporary=path.with_suffix('.tmp');tidy.to_parquet(temporary,index=False);temporary.replace(path)
        return tidy
    begin=code.index('    for attempt in range(1, MAX_ATTEMPTS + 1):')
    end=code.index('    gene_ids = tidy[',begin)
    code=code[:begin]+'    tidy = cached_tidy(model, interval, variant, scorer, variant_id)\n'+code[end:]
    sys.argv=[str(source/'code/make_preds.py'),str(target),str(mapping),model]
    exec(compile(code,str(source/'code/make_preds.py'),'exec'),{'__name__':'__main__','cached_tidy':cached_tidy})
    (out/'complete.json').write_text(json.dumps({'status':'complete','model':model,'variants':inputs.variant.nunique(),'pairs':len(inputs)},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variants',required=True);p.add_argument('--mapping',required=True)
    p.add_argument('--model',choices=['ALL_FOLDS','FOLD_1'],required=True)
    p.add_argument('--out',required=True);p.add_argument('--source',default='data/eqtl/source')
    p.add_argument('--cache',default='results/eqtl/fresh/cache');p.add_argument('--limit',type=int)
    run(**vars(p.parse_args()))
