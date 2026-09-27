"""Collect the streamed-prefilter large-database reruns into one JSON file.

Reads each run's run-summary.json and the `/usr/bin/time -v` block that the
rerun script writes to stderr.log. A multiple rerun with the release binary is
preferred; when that run is missing, the earlier streaming build's run is used
and its build is recorded, so the paper can say which binary each row came from.

    python3 analysis/scripts/collect_large_db_streamed.py
"""
import hashlib
import json
import re
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
ROOT = Path('/mnt/data1/sage-plus-scientific/large-db-20260925')
OUT = PAPER / 'analysis/data/large-db-streamed.json'
SOURCES = [
    ('runs-b13', 'v0.1.0-beta.13', ROOT / '../../build-cache/b13-rc/sage'),
    ('runs-v2-mp4', 'v0.1.0-beta.11 streaming build, executable SHA-256 prefix 47127fd7', ROOT / 'sage-streaming-v2'),
]
MULTIPLES = ['10x', '30x', '100x']


def wall_seconds(text):
    match = re.search(r'Elapsed \(wall clock\).*: ([\d:.]+)', text)
    parts = [float(p) for p in match.group(1).split(':')]
    seconds = 0.0
    for part in parts:
        seconds = seconds * 60 + part
    return seconds


def peak_gib(text):
    match = re.search(r'Maximum resident set size \(kbytes\): (\d+)', text)
    return int(match.group(1)) / 1048576


def sha256(path):
    path = Path(path).resolve()
    if not path.is_file():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect(multiple):
    name = f'hek-igc-{multiple}-prefilter'
    for folder, build, binary in SOURCES:
        run = ROOT / folder / name
        summary = run / 'run-summary.json'
        if not summary.is_file():
            continue
        data = json.loads(summary.read_text())
        stderr = (run / 'stderr.log').read_text()
        config = json.loads((run / 'config.json').read_text())
        return {
            'multiple': multiple,
            'name': name,
            'build': build,
            'binary_sha256': sha256(binary),
            'run_dir': str(run),
            'prefilter_min_matched_peaks': config['database'].get('prefilter_min_matched_peaks'),
            'max_memory_gb': config.get('max_memory_gb'),
            'wall_seconds': round(wall_seconds(stderr), 2),
            'peak_rss_gib': round(peak_gib(stderr), 3),
            'peptides_in_database': data['peptides_in_database'],
            'psms_at_one_percent_fdr': data['psms_at_one_percent_fdr'],
            'peptides_at_one_percent_fdr': data['peptides_at_one_percent_fdr'],
            'proteins_at_one_percent_fdr': data['proteins_at_one_percent_fdr'],
        }
    raise SystemExit(f'no completed run for {name}')


def main():
    rows = [collect(m) for m in MULTIPLES]
    OUT.write_text(json.dumps({'runs': rows}, indent=2) + '\n')
    for row in rows:
        print(row['multiple'], row['build'], row['wall_seconds'], row['peak_rss_gib'],
              row['psms_at_one_percent_fdr'])


if __name__ == '__main__':
    main()
