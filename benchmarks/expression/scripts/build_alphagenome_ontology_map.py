#!/usr/bin/env python3
"""Look up file_mappings in ENCODE and compare their ontology/assay with AG."""

from __future__ import annotations

import argparse
from collections import Counter
import subprocess
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import requests


ENCODE_BASE = "https://www.encodeproject.org"
SESSION = requests.Session()
ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "benchmarks/expression/data"
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_expression_mapping import match_samples


def encode_session():
    return SESSION


def load_cache(path: Path) -> dict[str, Any]:
    if path.exists():
        with path.open() as handle:
            return json.load(handle)
    return {}


def save_cache(path: Path, cache: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)
    tmp.replace(path)


def encode_get(accession: str, cache: dict[str, Any], retries: int, sleep: float) -> dict[str, Any]:
    if accession in cache:
        return cache[accession]

    if accession.startswith("ENCFF"):
        url = f"{ENCODE_BASE}/files/{accession}/?format=json"
    elif accession.startswith("ENCSR"):
        url = f"{ENCODE_BASE}/experiments/{accession}/?format=json"
    else:
        raise ValueError(f"Unsupported ENCODE accession: {accession}")

    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = encode_session().get(
                url,
                headers={"accept": "application/json"},
                timeout=(10, 60),
            )
            if response.status_code in (400, 401, 403, 404, 410, 422):
                raise ValueError(f"ENCODE {accession}: HTTP {response.status_code}")
            response.raise_for_status()
            cache[accession] = response.json()
            time.sleep(sleep)
            return cache[accession]
        except requests.RequestException as exc:
            last_error = exc
            time.sleep(sleep * attempt * 5)
    raise RuntimeError(f"Failed to fetch {accession}: {last_error}") from last_error


def get_nested(mapping: dict[str, Any], *keys: str) -> Any:
    value: Any = mapping
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def extract_ontology(payload: dict[str, Any]) -> tuple[str, str, str, str]:
    ontology = payload.get("biosample_ontology") or {}
    term_id = ontology.get("term_id") or ""
    term_name = ontology.get("term_name") or ""
    classification = ontology.get("classification") or ""
    assay = payload.get("assay_title") or payload.get("assay_term_name") or ""
    return term_id, term_name, classification, assay


def candidate_accessions(row: dict[str, str]) -> list[str]:
    candidates = [
        row.get("id", ""),
        row.get("targets_identifier_base", ""),
        row.get("original_id", ""),
    ]
    seen: set[str] = set()
    unique = []
    for accession in candidates:
        accession = accession.strip()
        if accession and accession not in seen and accession.startswith(("ENCFF", "ENCSR")):
            unique.append(accession)
            seen.add(accession)
    return unique


def read_sample_rows(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8-sig")
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("The sample list is empty")
    if all(line.startswith(("ENCFF", "ENCSR")) and not any(c in line for c in ",\t ") for line in lines):
        return [{"id": line} for line in dict.fromkeys(lines)]
    delimiter = "\t" if "\t" in lines[0] else ","
    rows = list(csv.DictReader(lines, delimiter=delimiter))
    if not rows or not any(key in rows[0] for key in ("id", "original_id", "targets_identifier_base")):
        raise ValueError("Use one ENCODE ID per line, or a CSV/TSV with an id column")
    return rows


def resolve_ontology(payload, cache, retries, sleep):
    """File records may keep biosample/assay metadata on their parent experiment."""
    own = extract_ontology(payload)
    if own[0] and own[3]:
        return own
    parent = payload.get("dataset") or get_nested(payload, "replicate", "experiment")
    if isinstance(parent, dict):
        term = extract_ontology(parent)
        if term[0] and term[3]:
            return term
        parent = parent.get("accession") or parent.get("@id", "")
    if isinstance(parent, str):
        import re
        match = re.search(r"ENCSR[A-Z0-9]+", parent)
        if match:
            return extract_ontology(encode_get(match.group(), cache, retries, sleep))
    return own


def read_file_mappings(path, organism):
    rows = read_sample_rows(Path(path))
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
    candidates=candidate_accessions(row)
    if not candidates:errors.append('No ENCFF/ENCSR accession found in id/original_id/targets_identifier_base')
    for accession in candidates:
        try:
            payload=encode_get(accession,cache,retries,sleep)
            term,name,classification,assay=resolve_ontology(payload,cache,retries,sleep)
            if term and assay:
                result.update(query_accession=accession,query_type='file' if accession.startswith('ENCFF') else 'experiment',
                              ontology_curie=term,biosample_name=name,biosample_classification=classification,assay=assay)
                break
            errors.append(f'{accession}: missing ontology or assay')
        except Exception as error:errors.append(f'{accession}: {error}')
    result['errors']=' | '.join(errors)
    return result

def resolve_rows(rows, cache_path, retries=3, sleep=0.1):
    cache_path = Path(cache_path)
    cache = load_cache(cache_path)
    cached_count = len(cache)
    result = []
    try:
        for index, row in enumerate(rows, 1):
            result.append(resolve_one((index, row), cache, retries, sleep))
            if index % 25 == 0:
                if len(cache) != cached_count:
                    save_cache(cache_path, cache)
                    cached_count = len(cache)
                print(f"ENCODE: {index}/{len(rows)}", flush=True)
    finally:
        if len(cache) != cached_count or not cache_path.exists():
            save_cache(cache_path, cache)
    return result

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
        subprocess.run([sys.executable,str(Path(__file__).with_name('dump_alphagenome_rna_metadata.py')),
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

def run(file_mappings, organism, output, cache, metadata_path=None, ontology_path=None):
    suffix = "_mouse" if organism == "mouse" else ""
    output = Path(output)
    ontology_path = Path(ontology_path) if ontology_path else output.parent / f"qnorm_alphagenome_ontology_map{suffix}.csv"
    metadata_path = Path(metadata_path) if metadata_path else DATA / f"alphagenome_{organism}_rna_seq_metadata.csv"
    rows = read_file_mappings(file_mappings, organism)
    resolved = resolve_rows(rows, cache)
    write_csv(ontology_path, resolved)
    metadata = get_ag_metadata(organism, metadata_path, reuse=metadata_path.exists())
    summary = support_rows(resolved, metadata)
    write_csv(output, summary)
    matched = [r for r in summary if r["alphagenome_ontology_assay_supported"] == "yes"]
    print(f"Saved {len(summary)} rows to {output}; {len(matched)} samples match "
          f"{len({r['alphagenome_ontology_name'] for r in matched})} AG ontology-assay names.")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file-mappings", "--samples", "--qnorm-map", dest="samples", required=True)
    parser.add_argument("--organism", choices=["human", "mouse"], default="human")
    parser.add_argument("--metadata", help="AG RNA-seq metadata CSV (default: the saved organism metadata)")
    parser.add_argument("--out", help="Resolved ENCODE ontology CSV")
    parser.add_argument("--support-summary", help="AG support CSV (default: benchmarks/expression/data/qnorm_alphagenome_support_summary[_mouse].csv)")
    parser.add_argument("--cache", help="Raw ENCODE response cache")
    args = parser.parse_args()
    suffix = "_mouse" if args.organism == "mouse" else ""
    run(args.samples, args.organism,
        args.support_summary or DATA / f"qnorm_alphagenome_support_summary{suffix}.csv",
        args.cache or ROOT / f"data/encode_metadata_cache{suffix}.json",
        args.metadata, args.out)


if __name__ == "__main__":
    main()
