#!/usr/bin/env python3
"""Match a user's ENCODE samples to exact AlphaGenome RNA-seq output names."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HELPERS = ROOT / 'benchmarks/expression/scripts'

def read_rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle))

def match_samples(samples, metadata):
    canonical_names = {}
    for row in metadata:
        ontology = row.get('ontology_curie', '').strip()
        name = row.get('name', '').strip()
        if name.startswith(ontology + ' ') and ontology:
            # AlphaGenome's display assay title can be broader than its output
            # name. Prefer the exact name used by score_interval and notebooks.
            canonical_names.setdefault((ontology, name[len(ontology):].strip().casefold()), set()).add(name)
    mappings = {}
    report = []
    for row in samples:
        sample = row.get('id', '').strip()
        term, assay = row.get('ontology_curie', '').strip(), row.get('assay', '').strip()
        key = (term, assay.casefold())
        candidates = canonical_names.get(key, set())
        if not sample:
            reason = 'missing_sample_id'
        elif not term:
            reason = 'unresolved_ontology'
        elif not assay:
            reason = 'unresolved_assay'
        elif not candidates:
            reason = 'no_matching_ontology_assay'
        elif len(candidates) != 1:
            reason = 'ambiguous_metadata_names'
        else:
            reason = 'matched'
            name = next(iter(candidates))
            group = mappings.setdefault(name, [])
            if sample not in group:
                group.append(sample)
        report.append({'id': sample, 'ontology_curie': term, 'biosample_name': row.get('biosample_name', ''),
                       'assay': assay, 'ontology_name': next(iter(candidates)) if len(candidates) == 1 else '',
                       'status': reason, 'errors': row.get('errors', '')})
    return mappings, report

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--samples', help='ENCFF/ENCSR IDs, one per line or CSV/TSV with id header')
    source.add_argument('--ontology-map', help='Reuse a resolved ENCODE ontology CSV; no ENCODE requests')
    parser.add_argument('--organism', choices=['human', 'mouse'], default='human')
    parser.add_argument('--metadata', help='Reuse an exported AlphaGenome RNA-seq metadata CSV; no AG request')
    parser.add_argument('--out', required=True, help='Directory for JSON mappings and matching report')
    args = parser.parse_args()
    out = Path(args.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    resolved = Path(args.ontology_map).resolve() if args.ontology_map else out / 'encode_ontology.csv'
    metadata = Path(args.metadata).resolve() if args.metadata else out / 'alphagenome_rna_metadata.csv'
    if not args.metadata:
        subprocess.run([sys.executable, str(HELPERS / 'dump_alphagenome_rna_metadata.py'),
                        '--organism', args.organism, '--out', str(metadata)], check=True)
    if args.samples:
        subprocess.run([sys.executable, str(HELPERS / 'build_alphagenome_ontology_map.py'),
                        '--samples', args.samples, '--organism', args.organism, '--out', str(resolved),
                        '--support-summary', str(out / 'support_summary.csv'), '--cache', str(out / 'encode_cache.json'),
                        '--metadata', str(metadata)], check=True)
    groups, report = match_samples(read_rows(resolved), read_rows(metadata))
    if not report:
        parser.exit(2, 'No input samples were found.\n')
    with (out / 'sample_matching.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report[0])); writer.writeheader(); writer.writerows(report)
    # A representative ID is only a label for the prediction. It is not GT aggregation.
    first = {name: ids[0] for name, ids in sorted(groups.items())}
    all_ids = {name: ids for name, ids in sorted(groups.items())}
    for name, value in [('track_to_first_id.json', first), ('track_to_all_ids.json', all_ids)]:
        (out / name).write_text(json.dumps(value, indent=2) + '\n')
    manifest = {'organism': args.organism, 'resolved_ontology_csv': str(resolved), 'metadata_csv': str(metadata),
                'sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [resolved, metadata]},
                'input_rows': len(report), 'matched_rows': sum(r['status'] == 'matched' for r in report),
                'ontology_assay_groups': len(groups), 'unique_ontology_ids': len({name.split()[0] for name in groups})}
    (out / 'mapping_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"Matched {manifest['matched_rows']}/{len(report)} samples to {len(groups)} ontology-assay names.")
    print(f"Review {out / 'sample_matching.csv'} for unmatched samples.")
    if not groups:
        parser.exit(2, 'No matching ontology-assay names; inspect the matching report.\n')

if __name__ == '__main__': main()
