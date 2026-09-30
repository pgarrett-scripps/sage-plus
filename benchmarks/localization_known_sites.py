#!/usr/bin/env python3
"""Measure PTM localization error against known phosphosites.

PXD000138 (Marx et al. 2013) is a synthetic phosphopeptide library in which
every peptide carries one known phosphorylation site. Each run searches
three of its files with phospho STY plus Met oxidation and writes every
localization (`localization_q_value` 1.0), once per `min_separating_margin`
value. A localization is evaluable when its peptide occurs in the library
with exactly one seeded site; it is correct when it names one site and that
site is the seed. The estimated FLR at a cutoff is the largest reported
localization q-value among the passing rows.

PXD007058 is searched once at the default margin to count phospho and
oxidation localizations and protein sites in a biological phosphoproteome.

usage: localization_known_sites.py --sage BIN --work DIR --output JSON
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
from pathlib import Path

import pandas as pd

from provenance import atomic_json, sha256

DATA = Path('/mnt/data1/sage-plus-scientific/loc-known-sites')
MARGINS = [0, 1, 2, 3]
CUTOFFS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]


def search(sage: Path, config: dict, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    path = output / 'config.json'
    path.write_text(json.dumps(config, indent=1))
    if not (output / 'run-summary.json').exists():
        environment = os.environ | {'RAYON_NUM_THREADS': '8'}
        with (output / 'search.log').open('w') as log:
            subprocess.run([str(sage), str(path), '--output_directory', str(output),
                            '--disable-telemetry-i-dont-want-to-improve-sage'],
                           check=True, stdout=log, stderr=subprocess.STDOUT, env=environment)
    return {'config_sha256': sha256(path),
            'sites_sha256': sha256(output / 'results.sage.ptm-sites.parquet')}


def localizations(run: Path) -> pd.DataFrame:
    sites = pd.read_parquet(run / 'results.sage.ptm-sites.parquet')
    grouped = sites.groupby(['psm_id', 'modification'], sort=False)
    return grouped.agg(peptide=('peptide', 'first'), q=('localization_q_value', 'first'),
                       k=('position', 'size'),
                       sites=('position', lambda x: frozenset(p - 1 for p in x))).reset_index()


def truth_lookup():
    library = json.loads((DATA / 'pxd000138/truth_libraries.json').read_text())
    cache: dict[str, tuple[bool, frozenset]] = {}

    def truth(sequence: str):
        if sequence not in cache:
            seeds, found = set(), False
            for entry in library.values():
                period, seed, full = entry['L'], entry['seed'], entry['seq']
                start = full.find(sequence)
                while start != -1:
                    found = True
                    first = start - start % period
                    for unit in range(first, start + len(sequence), period):
                        position = unit + seed - start
                        if 0 <= position < len(sequence):
                            seeds.add(position)
                    start = full.find(sequence, start + 1)
            cache[sequence] = (found, frozenset(seeds))
        return cache[sequence]
    return truth


def known_site_curve(run: Path, truth) -> dict:
    loc = localizations(run)
    loc = loc[loc.modification.str.contains('79.96')].copy()
    loc['sequence'] = loc.peptide.str.replace(r'\[[^\]]*\]', '', regex=True).str.replace('-', '')
    looked = loc.sequence.map(truth)
    loc['library'] = [found for found, _ in looked]
    loc['seed'] = [seeds for _, seeds in looked]
    evaluable = loc[loc.library & (loc.seed.map(len) == 1) & (loc.k == 1)].copy()
    evaluable['correct'] = [sites <= seed for sites, seed in zip(evaluable.sites, evaluable.seed)]
    points = []
    for cutoff in CUTOFFS:
        passing = evaluable[evaluable.q <= cutoff]
        wrong = int((~passing.correct).sum())
        points.append({'cutoff': cutoff, 'passing': len(passing), 'false': wrong,
                       'true_flr': wrong / max(len(passing), 1),
                       'estimated_flr': float(passing.q.max()) if len(passing) else None})
    return {'evaluable': len(evaluable), 'competing': int((evaluable.q < 1).sum()),
            'points': points}


def per_type(run: Path) -> dict:
    loc = localizations(run)
    protein = pd.read_parquet(run / 'results.sage.protein-sites.parquet')
    out = {}
    for modification in sorted(loc.modification.unique()):
        rows = loc[loc.modification == modification]
        sites = protein[protein.modification == modification]
        out[modification] = {
            'localizations': len(rows),
            'localizations_q01': int((rows.q <= 0.01).sum()),
            'certain': int((rows.q == 0).sum()),
            'protein_sites_q01': int(((sites.best_localization_q_value <= 0.01)
                                      & (sites.site_q_value <= 0.01)).sum()),
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sage', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    truth = truth_lookup()
    base = json.loads((DATA / 'pxd000138/config.json').read_text())
    runs, curves = {}, {}
    for margin in MARGINS:
        config = copy.deepcopy(base)
        config['ptm_localization']['min_separating_margin'] = margin
        run = args.work / f'pxd000138-margin{margin}'
        runs[run.name] = search(args.sage, config, run)
        curves[str(margin)] = known_site_curve(run, truth)
    biological = json.loads((DATA / 'pxd007058/config.json').read_text())
    run = args.work / 'pxd007058'
    runs[run.name] = search(args.sage, biological, run)
    example = pd.read_parquet(args.work / 'pxd000138-margin1/results.sage.parquet',
                              columns=['localized_peptide'])['localized_peptide'].dropna()
    atomic_json(args.output, {
        'sage_sha256': sha256(args.sage), 'margins': MARGINS, 'cutoffs': CUTOFFS,
        'known_sites': curves, 'pxd007058': per_type(run), 'runs': runs,
        'localized_peptide': {'rows': len(example),
                              'scored_groups': int(example.str.contains('#g').sum())},
        'truth': str(DATA / 'pxd000138/truth_libraries.json'),
    })


if __name__ == '__main__':
    main()
