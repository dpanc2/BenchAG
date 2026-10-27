#!/usr/bin/env python3
"""Run human and mouse expression predictions sequentially; safe for nohup."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-dir', required=True)
    args = parser.parse_args()
    env = os.environ.copy()
    # A user-provided dotenv file is read only for this key and never logged.
    key_file = Path(env.get('BENCHAG_API_ENV', str(ROOT / '.env'))).expanduser()
    if not env.get('ALPHAGENOME_API_KEY') and key_file.is_file():
        for line in key_file.read_text().splitlines():
            line = line.strip()
            if line.startswith('export '):
                line = line[7:].strip()
            name, sep, value = line.partition('=')
            if sep and name.strip() == 'ALPHAGENOME_API_KEY':
                env['ALPHAGENOME_API_KEY'] = value.strip().strip('\"\'')
    if not env.get('ALPHAGENOME_API_KEY'):
        parser.exit(2, 'No API key: set ALPHAGENOME_API_KEY or BENCHAG_API_ENV (default .env).\n')
    for name in ('HTTP_PROXY', 'HTTPS_PROXY', 'ALL_PROXY', 'http_proxy', 'https_proxy', 'all_proxy'):
        env.pop(name, None)
    env['PYTHONUNBUFFERED'] = '1'
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    import fcntl
    lock = (ROOT / 'results/expression-128kb.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.exit(2, 'A 128 kb driver is already running.\n')
    record = {'pid': os.getpid(), 'started_utc': datetime.now(timezone.utc).isoformat(),
              'sequence_length': 131072, 'gene_mask_width': 100001, 'model_version': 'FOLD_1',
              'splits': ['test', 'valid'], 'status': 'running', 'jobs': []}
    manifest = run_dir / 'run.json'
    def save():
        temporary = manifest.with_suffix('.tmp')
        temporary.write_text(json.dumps(record, indent=2) + '\n')
        temporary.replace(manifest)
    save()
    for organism in ('human', 'mouse'):
        config = ROOT / f'configs/expression.{organism}.128kb.yaml'
        log = run_dir / f'{organism}.log'
        cmd = [sys.executable, '-u', str(ROOT / 'scripts/run_expression.py'), '--config', str(config)]
        job = {'organism': organism, 'config': str(config), 'command': cmd, 'log': str(log), 'status': 'running'}
        record['jobs'].append(job)
        save()
        print(f'Starting {organism}: {log}', flush=True)
        with log.open('a') as stream:
            proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
            job['pid'] = proc.pid
            save()
            code = proc.wait()
        job.update(exit_code=code, status='completed' if code == 0 else 'failed', finished_utc=datetime.now(timezone.utc).isoformat())
        save()
        print(f'{organism} exit code: {code}', flush=True)
    record.update(status='completed' if all(j['exit_code'] == 0 for j in record['jobs']) else 'failed', finished_utc=datetime.now(timezone.utc).isoformat())
    save()
    return 0 if record['status'] == 'completed' else 1

if __name__ == '__main__':
    sys.exit(main())
