#!/usr/bin/env python3
"""Download and verify public ATAC inputs, reference and ChromBPNet model."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/atac'
BENCH=ROOT/'benchmarks/atac'
SOURCE_COMMIT='9e6a91157702b80efbc3b43f23c06e21f7bae856'
HELPER_COMMIT='09938fdb4397ec0006510e5251e48920a505d4de'
MD5={'ENCFF650PHA.bigWig':'6254c06f028aa67f99958c289c8adbed',
     'ENCFF113GSV.tar.gz':'f1dd105bcfa444e30b645464f4206fda',
     'ENCFF802KJC.tar.gz':'40bbb77a35f5b35cc5fd13fbfa78ee7f'}
REFERENCE={
    'chr1':{'length':248956422,'canonical_sha256':'209d95937d1fc25d970390c56f11f79ac9d5ac127ba48c21fbd4e1aac3bb60a5'},
    'chr3':{'length':198295559,'canonical_sha256':'1500ed6a7516eaa5645c56168ff7449f7aadfbd341b94992f63ccdf91b96977c'},
    'chr6':{'length':170805979,'canonical_sha256':'199d76944a1b320cc8e9b3935acd86df00e3bfc3dfa776e6f391f3fc65bf9c6f'}}


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
                # The source one-hot encoder maps every non-ACGT letter to zeros.
                bases=line.strip().upper().translate(bytes.maketrans(b'RYSWKMBDHV', b'NNNNNNNNNN'))
                seen[current][0].update(bases);seen[current][1]+=len(bases)
    actual={c:{'length':n,'canonical_sha256':h.hexdigest()} for c,(h,n) in seen.items()}
    if actual!=REFERENCE:raise ValueError('Reference model-input bases differ from the verified hg38 reference')
    return actual


def build_peak_bed(peaks, version, output):
    """Colab method: exclude whole peaks overlapping 1 Mb AG training intervals."""
    from alphagenome.data import fold_intervals
    from alphagenome.models import dna_client
    train = fold_intervals.get_fold_intervals(
        model_version=dna_client.ModelVersion[version],
        organism=dna_client.Organism.HOMO_SAPIENS,
        subset=fold_intervals.Subset.TRAIN,
        example_regions_path=str(BENCH/'data/sequences_human.bed.gz'))
    center = (train.start + train.end)//2
    train['start'], train['end'] = center-2**19, center+2**19
    overlap = np.zeros(len(peaks), dtype=bool)
    for chrom, group in peaks.groupby(0):
        intervals = train[train.chromosome==chrom].sort_values('start')
        if intervals.empty: continue
        starts = intervals.start.to_numpy()
        ends = np.maximum.accumulate(intervals.end.to_numpy())
        previous = np.searchsorted(starts,group[2].to_numpy(),side='left')-1
        overlap[group.index] = (previous>=0) & (ends[np.maximum(previous,0)]>group[1].to_numpy())
    filtered = peaks.loc[~overlap].sort_values([0,1])
    filtered.to_csv(output,sep='\t',header=False,index=False)
    return filtered


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    bed=BENCH/'data/fold1_test_regions.bed'
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
    peak_counts = {}
    for version, expected in [('FOLD_0',11860),('FOLD_1',13812)]:
        generated = build_peak_bed(original,version,DATA/f'{version.lower()}_test_regions.bed')
        assert len(generated)==expected, f'{version}: unexpected peak count'
        peak_counts[version] = len(generated)
        if version=='FOLD_1':
            pd.testing.assert_frame_equal(
                generated.sort_values([0,1,2,3]).reset_index(drop=True),
                selected.sort_values([0,1,2,3]).reset_index(drop=True),
                check_exact=False,rtol=1e-12,atol=1e-12)
    print(f'Peak BEDs rebuilt from {len(original)} original peaks: {peak_counts}',flush=True)
    helper=BENCH/'chrombpnet/reformat_chrombpnet_h5.py'
    python=DATA/'chrombpnet_env/bin/python'
    if not (models/'chrombpnet_recompiled.h5').exists():
        subprocess.run([str(python),str(helper),'-cnb',str(models/'model.chrombpnet_nobias.fold_0.ENCSR200OML.h5'),'-bm',str(models/'model.bias_scaled.fold_0.ENCSR200OML.h5'),'-o',str(models)],check=True)
    ref=reference(DATA/'hg38_chr1_chr3_chr6.fa')
    from alphagenome.data import fold_intervals
    from alphagenome.models import dna_client
    audit={}
    for model in [dna_client.ModelVersion.FOLD_0,dna_client.ModelVersion.FOLD_1]:
        train=fold_intervals.get_fold_names(model,fold_intervals.Subset.TRAIN)
        partitions=BENCH/'data/sequences_human.bed.gz'
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
              'local_code_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [BENCH/'alpha/alpha.py',BENCH/'chrombpnet/chrom_bw.py',BENCH/'chrombpnet/chrom_metrics.py',helper]},
              'reference_urls':{c:f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{c}.fa.gz' for c in REFERENCE},
              'training_window_overlaps':audit,
              'peak_filter':{'method':'Exclude peaks overlapping AG training intervals resized to 1048576 bp',
                             'colab':'https://colab.research.google.com/drive/1qjR2W1LpRbrb0zmGescaiUGe7Ft1k0kN',
                             'original_peaks':len(original),'retained_peaks':peak_counts,'FOLD_1_matches_colleague':True},
              'extracted_files':{p.name:digest(p) for p in models.glob('*') if p.is_file()}}
    (DATA/'preparation.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS: ATAC inputs, reference, original test peaks and model folds verified.',flush=True)

if __name__=='__main__':main()
