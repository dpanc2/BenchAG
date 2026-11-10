#!/usr/bin/env python3
"""Download and verify public ATAC inputs, reference and ChromBPNet model."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/atac'
SOURCE_COMMIT='9e6a91157702b80efbc3b43f23c06e21f7bae856'
HELPER_COMMIT='09938fdb4397ec0006510e5251e48920a505d4de'
MD5={'ENCFF650PHA.bigWig':'6254c06f028aa67f99958c289c8adbed',
     'ENCFF113GSV.tar.gz':'f1dd105bcfa444e30b645464f4206fda',
     'ENCFF802KJC.tar.gz':'40bbb77a35f5b35cc5fd13fbfa78ee7f'}
REFERENCE={
    'chr1':{'length':248956422,'upper_sha256':'68a7731692bde93299f5994c938d6ba4eaa2890eaabe1564b52d895a6f9ad78f'},
    'chr3':{'length':198295559,'upper_sha256':'75fd8635511d01a22193b4f395915a24b78195782cfc716d11cd3eccdeb8974f'},
    'chr6':{'length':170805979,'upper_sha256':'9893618330289268c7a104789a596f24661e661e6d878f51e72ee7643f0187b4'}}


def digest(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def fetch(url,path,expected=None,algorithm='sha256'):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        print('Downloading',path.name,flush=True)
        temp=path.with_suffix(path.suffix+'.part')
        for attempt in range(3):
            try:
                with requests.get(url,stream=True,timeout=(30,180)) as r:
                    r.raise_for_status()
                    with temp.open('wb') as f:
                        for block in r.iter_content(1024*1024):f.write(block)
                temp.replace(path);break
            except requests.RequestException:
                if attempt==2:raise
    if expected and digest(path,algorithm)!=expected:
        raise ValueError(f'Checksum mismatch: {path}; remove it and retry')
    return path


def reference(destination):
    """Public UCSC hg38 chromosomes; verify exact bases used by the predictor."""
    destination=Path(destination)
    if not destination.exists():
        temp=destination.with_suffix('.tmp')
        with temp.open('wb') as output:
            for chrom in REFERENCE:
                url=f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz'
                path=fetch(url,DATA/'reference'/f'{chrom}.fa.gz')
                with gzip.open(path,'rb') as f:shutil.copyfileobj(f,output)
        temp.replace(destination)
    seen={};current=None
    with destination.open('rb') as f:
        for line in f:
            if line.startswith(b'>'):
                current=line[1:].split()[0].decode();seen[current]=[hashlib.sha256(),0]
            else:
                bases=line.strip().upper();seen[current][0].update(bases);seen[current][1]+=len(bases)
    actual={c:{'length':n,'upper_sha256':h.hexdigest()} for c,(h,n) in seen.items()}
    if actual!=REFERENCE:raise ValueError('Reference bases differ from the verified hg38 reference')
    return actual


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    source=DATA/'source'
    if not source.exists():subprocess.run(['git','clone','https://github.com/GalinaK413/atac_bench_example',str(source)],check=True)
    assert subprocess.check_output(['git','-C',str(source),'status','--porcelain','--untracked-files=no'],text=True).strip()=='', 'Source checkout has local changes'
    subprocess.run(['git','-C',str(source),'checkout','--detach',SOURCE_COMMIT],check=True)
    bed=source/'alpha/fold1_test_regions.bed'
    assert digest(bed)=='31985063199bfbcbfa4d43d19d92e305e2edc3774f73f1f8df729aba9d2d89b9'
    for name,md5 in MD5.items():
        accession=name.split('.')[0]
        fetch(f'https://www.encodeproject.org/files/{accession}/@@download/{name}',DATA/name,md5,'md5')
    models=DATA/'models';models.mkdir(exist_ok=True)
    for archive,names in [('ENCFF113GSV.tar.gz',[
        'model.chrombpnet_nobias.fold_0.ENCSR200OML.h5','model.bias_scaled.fold_0.ENCSR200OML.h5']),
        ('ENCFF802KJC.tar.gz',['peaks.testset.fold_0.ENCSR200OML.bed.gz'])]:
        with tarfile.open(DATA/archive) as tar:
            for name in names:
                target=models/name
                if not target.exists():
                    with tar.extractfile('./fold_0/'+name) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
    selected=pd.read_csv(bed,sep='\t',header=None)
    original=pd.read_csv(models/'peaks.testset.fold_0.ENCSR200OML.bed.gz',sep='\t',header=None)
    originals=set(original.iloc[:,:3].itertuples(index=False,name=None))
    assert all(row in originals for row in selected.iloc[:,:3].itertuples(index=False,name=None))
    helper=fetch(f'https://raw.githubusercontent.com/kundajelab/chrombpnet/{HELPER_COMMIT}/chrombpnet/helpers/postprocessing/reformat_chrombpnet_h5.py',DATA/'reformat_chrombpnet_h5.py',
                 '6d05dcae0e829789163d9885009d4cce73cbe8bdc2c0d39566a7913bacbd66c1')
    python=DATA/'chrombpnet_env/bin/python'
    if not (models/'chrombpnet_recompiled.h5').exists():
        subprocess.run([str(python),str(helper),'-cnb',str(models/'model.chrombpnet_nobias.fold_0.ENCSR200OML.h5'),'-bm',str(models/'model.bias_scaled.fold_0.ENCSR200OML.h5'),'-o',str(models)],check=True)
    ref=reference(DATA/'hg38_chr1_chr3_chr6.fa')
    from alphagenome.data import fold_intervals
    from alphagenome.models import dna_client
    audit={}
    for model in [dna_client.ModelVersion.FOLD_0,dna_client.ModelVersion.FOLD_1]:
        train=fold_intervals.get_fold_names(model,fold_intervals.Subset.TRAIN)
        partitions=fetch('https://raw.githubusercontent.com/calico/borzoi/5c9358222b5026abb733ed5fb84f3f6c77239b37/data/sequences_human.bed.gz',DATA/'sequences_human.bed.gz')
        folds=pd.read_csv(partitions,sep='\t',header=None,names=['chrom','start','end','fold']);folds=folds[folds.fold.isin(train)]
        overlaps=0
        for chrom,group in selected.groupby(0):
            intervals=folds[folds.chrom==chrom].sort_values('start')
            if intervals.empty:continue
            starts=intervals.start.to_numpy();ends=np.maximum.accumulate(intervals.end.to_numpy())
            centers=((group[1]+group[2])//2).to_numpy();i=np.searchsorted(starts,centers+500,side='left')-1
            overlaps+=int(((i>=0)&(ends[np.maximum(i,0)]>centers-500)).sum())
        audit[model.name]=overlaps
    assert len(selected)==13812 and audit=={'FOLD_0':13550,'FOLD_1':0}
    manifest={'source_commit':SOURCE_COMMIT,'helper_commit':HELPER_COMMIT,'regions':len(selected),
              'input_md5':MD5,'reference':ref,
              'reference_urls':{c:f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{c}.fa.gz' for c in REFERENCE},
              'training_window_overlaps':audit,
              'extracted_files':{p.name:digest(p) for p in models.glob('*') if p.is_file()}}
    (DATA/'preparation.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: ATAC inputs, reference, original test peaks and model folds verified.',flush=True)

if __name__=='__main__':main()
