"""Derive the manuscript from completed release searches and preserved failures."""
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

from provenance import atomic_json, file_identities, sha256
from refresh_paper import BASE, ROOT, EXT, PLUS, UPSTREAM, read
import run_scientific as runner
from scientific_entrapment import evaluate

REPO = Path(__file__).resolve().parents[1]
PAPER = REPO / 'paper'
SELECTED = BASE / 'evidence-selected'
EXT_SELECTED = BASE / 'extension-selected'
RETRIES = BASE / 'retries'


def link(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        target.symlink_to(source.resolve(), target_is_directory=source.is_dir())


def select_matrix(source, target, retry_missing):
    target.mkdir(parents=True, exist_ok=True)
    completed = []
    failed = []
    for item in source.iterdir():
        if not item.is_dir() or not (item / 'result.json').exists():
            if item.name not in ('matrix-status.json',):
                link(item, target / item.name)
            continue
        result = read(item / 'result.json')
        chosen = item
        if result['status'] != 'complete' and retry_missing:
            job = copy.deepcopy(result['signature']['job'])
            retry_root = RETRIES / source.name
            print('Retrying missing analytical cell', source.name, job['id'], flush=True)
            repeated = runner.run_job(job, retry_root, 900)
            if repeated['status'] == 'complete':
                chosen = retry_root / job['id']
                result = repeated
                if source.name.startswith('entrapment-'):
                    table = chosen / ('results.sage.parquet' if job['engine'] == 'plus' else 'results.sage.tsv')
                    evaluate(table, source / 'paired.txt', chosen)
        link(chosen, target / item.name)
        if result['status'] == 'complete':
            completed.append(item.name)
        else:
            failed.append({'id': item.name, 'status': result['status']})
    if (source / 'plan.json').exists():
        atomic_json(target / 'matrix-status.json', {
            'expected_jobs': [j['id'] for j in read(source / 'plan.json')['jobs']],
            'completed': completed, 'failures': failed, 'pending': []})


def command(*args):
    print('Running', *map(str, args), flush=True)
    if PAPER in Path(args[0]).parents:
        environment = os.environ.copy()
        environment['UV_CACHE_DIR'] = '/tmp/sage-paper-uvcache'
        subprocess.run(['uv', 'run', '--offline', '--no-sync', *map(str, args)],
                       check=True, cwd=PAPER / 'analysis', env=environment)
    else:
        subprocess.run([sys.executable, *map(str, args)], check=True, cwd=REPO)


def main():
    assert read(BASE / 'refresh-status.json')['status'] == 'searches_finished'
    for item in ROOT.iterdir():
        if item.name != 'runs':
            link(item, SELECTED / item.name)
    for suite in sorted((ROOT / 'runs').iterdir()):
        retry_missing = suite.name.startswith(('entrapment-', 'public-comparison')) or suite.name == 'ptm'
        select_matrix(suite, SELECTED / 'runs' / suite.name, retry_missing)
    select_matrix(EXT, EXT_SELECTED, False)
    for name in ('lfq-upstream', 'lfq-plus'):
        if read(EXT / name / 'result.json')['status'] != 'complete':
            raise RuntimeError('LFQ requires a recorded retry before collecting results')
    offset_root = BASE / 'mass-offset'
    if not (offset_root / 'matrix.json').exists():
        command(REPO / 'benchmarks/run_mass_offset.py', '--root', offset_root,
                '--sage', PLUS, '--paired-fasta', '/data/sage-plus-scientific/mass-offset-20260919/inputs/paired.fasta')
    command(REPO / 'benchmarks/run_named_modifications.py', '--sage', PLUS,
            '--output', PAPER / 'analysis/data/named-modifications.json', '--work', BASE / 'named-modifications')
    output = REPO / 'benchmarks/scientific-results/20260920'
    output.mkdir(parents=True, exist_ok=True)
    command(REPO / 'benchmarks/analyze_scientific.py', '--root', SELECTED,
            '--output', output / 'pilot-summary.json')
    offset_summary = REPO / 'benchmarks/scientific-results/mass-offset-20260920/summary.json'
    command(REPO / 'benchmarks/summarize_mass_offset.py', '--root', offset_root,
            '--pairs', '/data/sage-plus-scientific/mass-offset-20260919/inputs/paired.txt',
            '--output', offset_summary)
    context = {'evidence': str(SELECTED), 'extension': str(EXT_SELECTED),
               'analysis_repository': str(REPO),
               'sage_version': 'v0.15.0-beta.2', 'sage_plus_version': 'v0.1.0-beta.6',
               'sage_plus_commit': '3e30135fb8786ec8a12c1f62e0ff9300e57f9567',
               'executables': file_identities([UPSTREAM, PLUS])}
    atomic_json(PAPER / 'analysis/data/release-context.json', context)
    atomic_json(PAPER / 'analysis/data/report-extension/plan.json', read(EXT_SELECTED / 'plan.json'))
    for section in ('public', 'ptm', 'lfq', 'scaling'):
        command(PAPER / 'analysis/scripts/collect_report.py', section)
    command(PAPER / 'analysis/scripts/collect_matched_fdp.py')
    command(PAPER / 'analysis/scripts/audit_report.py')
    attempts = []
    paths = list((ROOT / 'runs').glob('*/*/result.json')) + list(EXT.glob('*/result.json')) + list(RETRIES.glob('*/*/result.json'))
    for path in sorted(paths):
        record = read(path)
        job = record['signature']['job']
        suite = job.get('suite', path.parent.parent.name)
        if suite.startswith('entrapment-'):
            suite = 'entrapment'
        attempts.append({'suite': suite, 'engine': job['engine'], 'id': job['id'],
                         'status': record['status'], 'record': str(path),
                         'retry': RETRIES in path.parents, 'exit_code': record.get('exit_code')})
    atomic_json(PAPER / 'analysis/data/execution.json', {'attempts': attempts,
                'inputs_sha256': file_identities(paths),
                'policy': 'One retry is allowed for a missing identification or entrapment cell. Failed attempts remain in the raw record. Timing and scaling do not replace failed trials.'})
    checked = {}
    for path in paths:
        record = read(path)
        for name, expected in (record['signature']['inputs'] | record.get('outputs', {})).items():
            if name not in checked:
                checked[name] = sha256(Path(name))
            assert checked[name] == expected, f'Changed evidence {name}'
    atomic_json(output / 'evidence-verification.json', {
        'status': 'verified', 'members_verified': len(checked), 'inputs_sha256': checked,
        'scope': 'All comparative attempt inputs and recorded outputs, including failed attempts.'})
    atomic_json(BASE / 'analysis-status.json', {'status': 'ready_for_manuscript_build',
                'attempts': len(attempts), 'failures': [r for r in attempts if r['status'] != 'complete']})


if __name__ == '__main__':
    main()
