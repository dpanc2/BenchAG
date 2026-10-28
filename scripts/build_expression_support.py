#!/usr/bin/env python3
"""file_mappings CSV -> ENCODE metadata -> AlphaGenome support-summary table."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HELPERS = ROOT / 'benchmarks/expression/scripts'
spec = importlib.util.spec_from_file_location('encode_support_resolver', HELPERS/'build_alphagenome_ontology_map.py')
encode = importlib.util.module_from_spec(spec);spec.loader.exec_module(encode)
# Reuse the exact ontology-and-assay matching used by the public prediction workflow.
sys.path.insert(0,str(ROOT/'scripts'))
from prepare_expression_mapping import match_samples

def read_file_mappings(path, organism):
    rows = encode.read_sample_rows(Path(path))
    columns = set(rows[0])
    column = next((c for c in ('genome','assembly') if c in columns), None)
    if column:
        builds = {'human':{'hg19','hg38','grch37','grch38'},
                  'mouse':{'mm9','mm10','mm39','grcm37','grcm38','grcm39'}}[organism]
        rows = [r for r in rows if r.get(column,'').strip().casefold() in builds]
        if not rows:raise ValueError(f'No {organism} rows found in {column}; check the genome labels')
    ids = [r.get('id','').strip() for r in rows]
    if any(not i for i in ids):raise ValueError('Every file-mapping row needs an id')
    if len(set(ids)) != len(ids):raise ValueError('Duplicate sample IDs in file mappings; resolve them first')
    return rows

def resolve_one(item, cache, retries, sleep):
    index,row = item
    result={'row_index':str(index),'id':row.get('id',''),'original_id':row.get('original_id',''),
            'targets_identifier_base':row.get('targets_identifier_base',''),
            'query_accession':'','query_type':'','ontology_curie':'','biosample_name':'',
            'biosample_classification':'','assay':'','errors':''}
    errors=[]
    candidates=encode.candidate_accessions(row)
    if not candidates:errors.append('No ENCFF/ENCSR accession found in id/original_id/targets_identifier_base')
    for accession in candidates:
        try:
            payload=encode.encode_get(accession,cache,retries,sleep)
            term,name,classification,assay=encode.resolve_ontology(payload,cache,retries,sleep)
            if term and assay:
                result.update(query_accession=accession,query_type='file' if accession.startswith('ENCFF') else 'experiment',
                              ontology_curie=term,biosample_name=name,biosample_classification=classification,assay=assay)
                break
            errors.append(f'{accession}: missing ontology or assay')
        except Exception as error:errors.append(f'{accession}: {error}')
    result['errors']=' | '.join(errors)
    return result

def resolve_rows(rows, cache_path, workers=2, retries=3, sleep=0.1):
    cache_path=Path(cache_path);cache=encode.load_cache(cache_path)
    output=[]
    cached_count=len(cache)
    by_index={}
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            tasks={pool.submit(resolve_one,item,cache,retries,sleep):item[0] for item in enumerate(rows,1)}
            for future in as_completed(tasks):
                result=future.result()
                by_index[tasks[future]]=result
                output.append(result)
                if len(output)%25==0:
                    if len(cache) != cached_count:
                        encode.save_cache(cache_path,cache.copy())
                        cached_count=len(cache)
                    print(f'ENCODE: {len(output)}/{len(rows)} rows processed',flush=True)
    finally:
        if len(cache) != cached_count or not cache_path.exists():
            encode.save_cache(cache_path,cache.copy())
    return [by_index[i] for i in sorted(by_index)]

def read_csv(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream))

def write_csv(path, rows):
    if not rows:raise ValueError('Cannot write an empty support table')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
    temporary.replace(path)

def get_ag_metadata(organism, path, reuse=False):
    path=Path(path)
    if not reuse:
        subprocess.run([sys.executable,str(HELPERS/'dump_alphagenome_rna_metadata.py'),
                        '--organism',organism,'--out',str(path)],check=True)
    return read_csv(path)

def support_rows(resolved, metadata):
    if not metadata or not {'name','ontology_curie','strand'}<=set(metadata[0]):
        raise ValueError('AG metadata needs name, ontology_curie and strand columns')
    _,matches=match_samples(resolved,metadata)
    ontology_counts=Counter(r['ontology_curie'] for r in metadata)
    name_counts=Counter(r['name'] for r in metadata)
    result=[]
    for row,match in zip(resolved,matches):
        term=row['ontology_curie'];name=match['ontology_name'] if match['status']=='matched' else ''
        ontology_support=('yes' if ontology_counts[term] else 'no') if term else 'unknown'
        result.append({**row,
            'alphagenome_rna_track_count':ontology_counts[term] if term else 0,
            'alphagenome_supported':ontology_support,
            'alphagenome_ontology_assay_supported':'yes' if name else ('unknown' if not term or not row['assay'] else 'no'),
            'alphagenome_ontology_name':name,
            'alphagenome_ontology_assay_track_count':name_counts[name] if name else 0,
            'alphagenome_strands':','.join(sorted({r['strand'] for r in metadata if name and r['name']==name})),
            'match_status':match['status']})
    return result

def run(file_mappings, organism, output, cache, metadata_path=None, workers=2):
    output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    selected=read_file_mappings(file_mappings,organism)
    print(f'Input: {len(selected)} {organism} sample rows',flush=True)
    resolved=resolve_rows(selected,cache,workers=workers)
    suffix='_mouse' if organism=='mouse' else ''
    ontology_path=output.parent/f'qnorm_alphagenome_ontology_map{suffix}.csv'
    write_csv(ontology_path,resolved)
    metadata_path=Path(metadata_path) if metadata_path else output.parent/f'alphagenome_{organism}_rna_seq_metadata.csv'
    metadata=get_ag_metadata(organism,metadata_path,reuse=metadata_path.exists())
    rows=support_rows(resolved,metadata);write_csv(output,rows)
    manifest={'organism':organism,'created_utc':datetime.now(timezone.utc).isoformat(),
        'file_mappings':str(Path(file_mappings).resolve()),'encode_cache':str(Path(cache).resolve()),
        'resolved_ontology_table':str(ontology_path),'ag_metadata':str(metadata_path.resolve()),'output':str(output),
        'selected_rows':len(selected),'ontology_resolved':sum(bool(r['ontology_curie']) for r in rows),
        'ontology_supported':sum(r['alphagenome_supported']=='yes' for r in rows),
        'exact_name_supported':sum(r['alphagenome_ontology_assay_supported']=='yes' for r in rows),
        'shared_names':len({r['alphagenome_ontology_name'] for r in rows if r['alphagenome_ontology_name']}),
        'sha256':{str(Path(p).resolve()):hashlib.sha256(Path(p).read_bytes()).hexdigest()
                  for p in [file_mappings,cache,metadata_path,ontology_path,output]}}
    output.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f"Saved {len(rows)} rows to {output}; {manifest['exact_name_supported']} samples match {manifest['shared_names']} AG ontology-assay names.",flush=True)
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--file-mappings',required=True,help='Original CSV/TSV with id and optional original_id, targets_identifier_base, genome')
    p.add_argument('--organism',choices=['human','mouse'],default='human')
    p.add_argument('--out',help='Default: benchmarks/expression/data/qnorm_alphagenome_support_summary[_mouse].csv')
    p.add_argument('--cache',help='Raw ENCODE API-response cache; missing entries are fetched from ENCODE')
    p.add_argument('--metadata',help='Existing AG RNA-seq metadata CSV; if omitted, export when the default file is absent')
    p.add_argument('--workers',type=int,choices=[1,2,3,4],default=2)
    a=p.parse_args();suffix='_mouse' if a.organism=='mouse' else ''
    output=Path(a.out) if a.out else ROOT/f'benchmarks/expression/data/qnorm_alphagenome_support_summary{suffix}.csv'
    cache=Path(a.cache) if a.cache else ROOT/f'data/encode_metadata_cache{suffix}.json'
    run(a.file_mappings,a.organism,output,cache,a.metadata,a.workers)

if __name__=='__main__':
    try:main()
    except (ValueError,FileNotFoundError) as error:raise SystemExit(f'Input error: {error}')
