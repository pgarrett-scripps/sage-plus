"""Refresh the manuscript comparison using the published Sage Plus executable."""
import copy
import json
import platform
import sys
from pathlib import Path

from provenance import atomic_json, sha256
import run_scientific as runner
from scientific_entrapment import evaluate

BASE = Path(__file__).resolve().parents[2]
OLD = Path('/data/sage-plus-scientific/20260914')
OLD_EXT = Path('/data/sage-plus-scientific/report-extension-20260915')
ROOT = BASE / 'evidence'
EXT = BASE / 'extension'
PLUS = BASE / 'bin/sage-plus-v0.1.0-beta.6-x86_64-unknown-linux-gnu/sage'
UPSTREAM = Path('/home/ty/Repos/sage-plus/benchmarks/.work/targets/baseline-df9219951cc9a54c/release/sage')
runner.BASELINE = UPSTREAM
runner.CANDIDATE = PLUS


def read(path):
    return json.loads(path.read_text())


def matrix(name, jobs, destination):
    destination.mkdir(parents=True, exist_ok=True)
    jobs = copy.deepcopy(jobs)
    for job in jobs:
        job['binary'] = str(PLUS if job['engine'] == 'plus' else UPSTREAM)
    plan = {'comparison': {'upstream': 'v0.15.0-beta.2', 'plus': 'v0.1.0-beta.6'}, 'jobs': jobs}
    atomic_json(destination / 'plan.json', plan)
    atomic_json(destination / 'environment.json', {
        'platform': platform.platform(), 'python': platform.python_version(),
        'baseline_sha256': sha256(UPSTREAM), 'candidate_sha256': sha256(PLUS),
        'cpu': Path('/proc/cpuinfo').read_text(), 'timeout_seconds_per_job': 900})
    completed = []
    failures = []
    for index, job in enumerate(jobs):
        print(name, index + 1, '/', len(jobs), job['id'], flush=True)
        result = runner.run_job(job, destination, 900)
        if result['status'] == 'complete':
            completed.append(job['id'])
        else:
            failures.append({'id': job['id'], 'status': result['status']})
        atomic_json(destination / 'matrix-status.json', {
            'expected_jobs': [j['id'] for j in jobs], 'completed': completed,
            'failures': failures, 'pending': [j['id'] for j in jobs[index + 1:]]})
    return failures


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    for name in ('references', 'ptm-truth', 'inputs', 'converted'):
        target = ROOT / name
        if not target.exists():
            target.symlink_to(OLD / name, target_is_directory=True)
    ptm = [j for j in read(OLD / 'runs/ptm/plan.json')['jobs'] if j['id'].endswith('-all')]
    matrix('ptm', ptm, ROOT / 'runs/ptm')
    matrix('public-timing', read(OLD / 'public-repeat-v2-plan.json')['jobs'], ROOT / 'runs/public-timing')
    matrix('extension', read(OLD_EXT / 'plan.json')['jobs'], EXT)
    public = [j for j in read(OLD / 'public-comparison-plan.json')['jobs'] if j['study'] == 'PXD001468']
    public += read(OLD / 'public-comparison-v2-plan.json')['jobs']
    matrix('public-comparison-v2', public, ROOT / 'runs/public-comparison-v2')
    for old in sorted((OLD / 'runs').glob('entrapment-*')):
        destination = ROOT / 'runs' / old.name
        destination.mkdir(parents=True, exist_ok=True)
        for name in ('paired.txt', 'paired.fasta', 'generation.json'):
            target = destination / name
            if not target.exists():
                target.symlink_to(old / name)
        jobs = [read(path)['signature']['job'] for path in sorted(old.glob('*/result.json'))]
        matrix(old.name, jobs, destination)
        for job in jobs:
            output = destination / job['id']
            if read(output / 'result.json')['status'] == 'complete' and not (output / 'calibration.json').exists():
                source = output / ('results.sage.parquet' if job['engine'] == 'plus' else 'results.sage.tsv')
                evaluate(source, destination / 'paired.txt', output)
    local = [j for j in read(OLD / 'runs/local-paired/plan.json')['jobs'] if j['workload'] != 'broad-ptm']
    matrix('local-paired', local, ROOT / 'runs/local-paired')
    atomic_json(BASE / 'refresh-status.json', {'status': 'searches_finished', 'root': str(ROOT), 'extension': str(EXT)})


if __name__ == '__main__':
    main()
