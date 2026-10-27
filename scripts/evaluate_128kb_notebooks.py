#!/usr/bin/env python3
"""Execute AlphaGenome-only comparison notebook cells against saved 128 kb CSVs."""
import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import sys

SOURCE_ROOT=Path('/home/jovyan/dpanc/GENA_LM/GENA_LM_expression_branch')
DATA_ROOT=Path('/home/jovyan/dpanc/benchmarking/data')
BENCHMARK_ROOT=Path('/home/jovyan/dpanc/benchmarking')

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    os.environ.update(GENA_HOME=str(SOURCE_ROOT),BENCHMARK_ROOT=str(BENCHMARK_ROOT),TASK_ROOT=str(BENCHMARK_ROOT),DATA_ROOT=str(DATA_ROOT))
    base=SOURCE_ROOT/'downstream_tasks/expression_prediction/benchmarks/alphagenome_benchmark'
    summary=[];manifests=[]
    for organism,runid in [('human','f01d5f895d9c85de'),('mouse','399f3ae01998ea32')]:
        for normalization in ['qnorm','without_qnorm']:
            name='gena_ag_gt_ontology_compare'+('_mm10' if organism=='mouse' else '')+('_no_qnorm' if normalization=='without_qnorm' else '')+'.ipynb'
            source=base/'notebooks'/name
            original=json.loads(source.read_text())
            target=out/f'{organism}_{normalization}';target.mkdir(exist_ok=True)
            shutil.copy2(source,target/'source.ipynb')
            truth=(DATA_ROOT/('human_id_full_tpms.csv' if organism=='human' else 'mouse_id_full_tpms.csv')) if normalization=='without_qnorm' else DATA_ROOT/('borzoi_and_v1_human_qnorm_matrix_by_id.csv' if organism=='human' else 'true_mouse_all_genes_qnorm_samples_by_genes.WITH_TEST.csv')
            mapping=base/'data'/('alphagenome_track_to_qnorm_ids_all.json' if organism=='human' else 'alphagenome_track_to_qnorm_ids_mouse_all.json')
            pred=out/'predictions'/runid
            setup=f'''import os, json, re
from pathlib import Path
import numpy as np
import pandas as pd
REPO_ROOT=Path({str(SOURCE_ROOT)!r})
BENCHMARK_ROOT=Path({str(BENCHMARK_ROOT)!r})
DATA_ROOT=Path({str(DATA_ROOT)!r})
TRUE_PATH=Path({str(truth)!r})
ONTOLOGY_JSON=Path({str(mapping)!r})
AG_DIR=Path({str(pred)!r})
OUTPUT_DIR=Path({str(target)!r})
MODEL_DIRS={{}}
SPLITS=['valid','test']
'''
            selected=[3,7] if organism=='human' else [5,6]
            codes=[setup]+[''.join(original['cells'][i]['source']) for i in selected]
            # The current human base notebook points at non-qnorm data by default;
            # TRUE_PATH is explicitly restored to the intended qnorm matrix above.
            # Aggregation, log transform and finite-pair Pearson helpers are unchanged.
            compute='''rows=[]
per_gene_rows=[]
per_track_rows=[]
for split in SPLITS:
    path=AG_DIR/f'predictions_{split}_intervals.csv'
    prediction=read_prediction(path)
    if prediction.gene_id.duplicated().any():raise ValueError('Duplicate prediction gene IDs')
    ag=GROUP_PREDICTION(prediction)
    truth=GROUPED_TRUTH
    if truth.index.duplicated().any():raise ValueError('Duplicate truth gene IDs')
    genes=truth.index.intersection(ag.index)
    all_groups=truth.columns.intersection(ag.columns)
    usable=[g for g in all_groups if np.isfinite(ag.loc[genes,g].to_numpy(dtype=float)).any()]
    if not len(genes) or not usable:raise ValueError('Empty aligned population')
    xdf=truth.loc[genes,usable];ydf=ag.loc[genes,usable]
    x=xdf.to_numpy(dtype=float);y=ydf.to_numpy(dtype=float)
    per_track=CORRELATE(x,y,axis=0)
    per_gene=CORRELATE(x,y,axis=1)
    gene_corr=float(np.nanmean(per_track))
    cell_corr=float(np.nanmean(per_gene))
    assert np.isfinite(gene_corr) and np.isfinite(cell_corr)
    np.testing.assert_allclose(mean_correlations(xdf,ydf),[gene_corr,cell_corr],rtol=1e-12)
    # Independently check representative finite-pair correlations.
    for axis,values in [(0,per_track),(1,per_gene)]:
        checked=0
        for i,r in enumerate(values):
            if not np.isfinite(r):continue
            xx=x[:,i] if axis==0 else x[i,:]
            yy=y[:,i] if axis==0 else y[i,:]
            finite=np.isfinite(xx)&np.isfinite(yy)
            np.testing.assert_allclose(r,np.corrcoef(xx[finite],yy[finite])[0,1],rtol=1e-8,atol=1e-10)
            checked+=1
            if checked==5:break
    available=(np.isfinite(x)&np.isfinite(y)).sum(axis=1)
    rows.append(dict(organism=ORGANISM,normalization=NORMALIZATION,split=split,sequence_length=131072,
        gene_corr=gene_corr,cell_corr=cell_corr,n_genes=len(genes),n_track_groups=len(usable),
        n_aligned_groups_including_unavailable=len(all_groups),n_unique_ontology_curies=len({str(g).split(maxsplit=1)[0] for g in usable}),
        finite_gene_correlations=int(np.isfinite(per_gene).sum()),finite_track_correlations=int(np.isfinite(per_track).sum()),
        mean_finite_groups_per_gene=float(available.mean()),min_finite_groups_per_gene=int(available.min()),max_finite_groups_per_gene=int(available.max())))
    per_gene_rows.extend(dict(split=split,gene_id=str(g),cell_corr=float(r)) for g,r in zip(genes,per_gene))
    per_track_rows.extend(dict(split=split,track_group=str(g),gene_corr=float(r)) for g,r in zip(usable,per_track))
pd.DataFrame(rows).to_csv(OUTPUT_DIR/'summary.csv',index=False,float_format='%.17g')
pd.DataFrame(per_gene_rows).to_csv(OUTPUT_DIR/'per_gene_correlations.csv',index=False,float_format='%.17g')
pd.DataFrame(per_track_rows).to_csv(OUTPUT_DIR/'per_track_correlations.csv',index=False,float_format='%.17g')
print(pd.DataFrame(rows).to_string(index=False))
'''
            aliases=('GROUPED_TRUTH=true_grouped\nGROUP_PREDICTION=lambda df: aggregate_qnorm_ids(df,TRACK_TO_QNORM_IDS)[0]\nCORRELATE=pairwise_corr\n' if organism=='human' else 'GROUPED_TRUTH=true_ontology\nGROUP_PREDICTION=lambda df: mean_by_ontology(df,ontology_to_ids)\nCORRELATE=_pairwise_corr\n')
            codes.append(f'ORGANISM={organism!r}\nNORMALIZATION={normalization!r}\n'+aliases+compute)
            notebook={'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'execution_method':'Inspected cells executed sequentially in one Python namespace; nbclient unavailable','source_notebook':str(source),'source_cell_indices':selected},'cells':[]}
            namespace={'__name__':'__main__'}
            for count,code in enumerate(codes,1):
                cell={'cell_type':'code','id':f'cell-{count}','metadata':{},'source':code.splitlines(keepends=True),'execution_count':count,'outputs':[]}
                buffer=io.StringIO()
                try:
                    with contextlib.redirect_stdout(buffer),contextlib.redirect_stderr(buffer):exec(compile(code,f'{name}:execution-cell-{count}','exec'),namespace)
                finally:
                    if buffer.getvalue():cell['outputs'].append({'output_type':'stream','name':'stdout','text':buffer.getvalue().splitlines(keepends=True)})
                    notebook['cells'].append(cell)
                    (target/'executed.ipynb').write_text(json.dumps(notebook,indent=1))
            summary.extend(namespace['rows'])
            manifest={'organism':organism,'normalization':normalization,'source_notebook':str(source),'source_sha256':digest(source),'truth':str(truth),'truth_sha256':digest(truth),'mapping':str(mapping),'mapping_sha256':digest(mapping),'predictions':{s:{'path':str(pred/f'predictions_{s}_intervals.csv'),'sha256':digest(pred/f'predictions_{s}_intervals.csv')} for s in ['valid','test']},'execution_method':notebook['metadata']['execution_method'],'source_cell_indices':selected,'changes':['128 kb prediction paths','AlphaGenome only; GENA loops/exports omitted','Explicit qnorm/non-qnorm truth paths','Unambiguous axis labels correcting reversed human export labels','Entirely missing track groups omitted; original finite-pair mask and n>3 rule unchanged'],'python':sys.executable,'finished_utc':datetime.now(timezone.utc).isoformat()}
            (target/'manifest.json').write_text(json.dumps(manifest,indent=2));manifests.append(manifest)
            print(organism,normalization,namespace['rows'],flush=True)
    import pandas as pd
    pd.DataFrame(summary).to_csv(out/'summary.csv',index=False,float_format='%.17g')
    (out/'manifest.json').write_text(json.dumps({'status':'completed','results':manifests,'metric_definitions':{'gene_corr':'Mean Pearson across genes, averaged over track groups','cell_corr':'Mean Pearson across track groups, averaged over genes'},'source_script_sha256':digest(__file__)},indent=2))
    print('COMPLETED',out,flush=True)
if __name__=='__main__':main()
