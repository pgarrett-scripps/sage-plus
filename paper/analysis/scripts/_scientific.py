"""Read the frozen Sage and Sage Plus release comparison."""
import json
from pathlib import Path
from statistics import median

PAPER = Path(__file__).resolve().parents[2]
PILOT = '../benchmarks/scientific-results/20260914/pilot-summary.json'
ARCHIVE = '../benchmarks/scientific-results/20260914/archive-verification.json'
INPUTS = [PILOT, ARCHIVE]
ENGINE = {'upstream': 'Sage', 'plus': 'Sage Plus'}
STUDY = {'human': 'HEK', 'hye': 'Mixture'}


def load():
    return [json.loads((PAPER / path).read_text()) for path in INPUTS]


def timing(pilot):
    rows = []
    for study in ('PXD001468', 'PXD028735'):
        for engine in ('upstream', 'plus'):
            group = [j for j in pilot['jobs'] if j['suite'] == 'public-timing'
                     and j['study'] == study and j['engine'] == engine and not j['warmup']]
            assert len(group) == 3 and all(j['status'] == 'complete' for j in group)
            rows.append(dict(study=study, engine=engine, trials=len(group),
                             seconds=median(j['wall_seconds'] for j in group),
                             seconds_min=min(j['wall_seconds'] for j in group),
                             seconds_max=max(j['wall_seconds'] for j in group),
                             rss=median(j['peak_rss_mib'] for j in group),
                             psms=median(j['target_psms'] for j in group)))
    return rows


def add_stats(st):
    pilot, archive = load()
    for key, value, desc in (
        ('archive.members', archive['members_verified'], 'Individually verified archive members'),
        ('pilot.bootstrap', pilot['calibration_uncertainty'][0]['bootstrap_replicates'], 'Paired bootstrap resamples'),
        ('pilot.public_pairs', sum(r['status'] == 'complete' for r in pilot['public_identifications']), 'Completed public identification pairs'),
    ):
        st.add(key, value, fmt=',', desc=desc, sign='+')
    timed = timing(pilot)
    for row in timed:
        for key, fmt in (('seconds', '.2f'), ('rss', '.1f')):
            st.add(f"pilot.{row['study']}.{row['engine']}.{key}", row[key], fmt=fmt,
                   desc=f"{row['study']} {ENGINE[row['engine']]} timing median {key}", sign='+')
    for study in ('PXD001468', 'PXD028735'):
        rows = {r['engine']: r for r in timed if r['study'] == study}
        for key in ('seconds', 'rss'):
            reduction = 100 * (1 - rows['plus'][key] / rows['upstream'][key])
            st.add(f'pilot.{study}.{key}_reduction', reduction, fmt='.1f',
                   desc=f'{study} Sage Plus reduction in {key} relative to Sage', sign='+', between=(0,100))
    for row in pilot['calibration_uncertainty']:
        for engine, value in row['means'].items():
            st.add(f"pilot.{row['study']}.{engine}.fdp", value * 100, fmt='.3f',
                   desc=f"{STUDY[row['study']]} {ENGINE[engine]} mean paired FDP percent", between=(0,100))
        st.add(f"pilot.{row['study']}.delta", (row['means']['plus'] - row['means']['upstream']) * 100,
               fmt='+.3f', desc='Sage Plus minus Sage FDP in percentage points', between=(-100,100))
        for bound, value in zip(('lower', 'upper'), row['percentile_intervals']['paired_difference']):
            st.add(f"pilot.{row['study']}.{bound}", value * 100, fmt='+.3f',
                   desc='Conditional difference interval in percentage points', between=(-100,100))
    overlap = [r['jaccard'] for r in pilot['public_identifications'] if r['status'] == 'complete']
    for name, value in (('min', min(overlap)), ('max', max(overlap))):
        st.add(f'pilot.overlap.{name}', value, fmt='.4f', desc=f'{name} public PSM Jaccard', between=(0,1))
