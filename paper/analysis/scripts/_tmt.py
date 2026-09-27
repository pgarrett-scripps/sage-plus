"""Read the HYpro16 TMT ratio-compression summary (PXD020815)."""
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
SUMMARY = 'analysis/data/tmt-hypro16.json'
TMT_INPUTS = [SUMMARY]
RUNS = {'upstream': 'upstream-{method}-mzml', 'plus': 'plus-b13-{method}-mzml',
        'plus.raw': 'plus-b13-{method}-raw'}


def load():
    return json.loads((PAPER / SUMMARY).read_text())


def add_tmt_stats(st):
    runs = load()['runs']
    for method in ('hrms2', 'sps'):
        for engine, pattern in RUNS.items():
            row = runs[pattern.format(method=method)]
            st.add(f'tmt.{method}.{engine}.fold20', row['levels']['1/20']['observed_fold_1x_over_level'],
                   fmt='.3f', desc=f'{engine} {method} observed yeast fold at the 20:1 level', sign='+')
            st.add(f'tmt.{method}.{engine}.slope', row['compression_slope_spiked_levels'],
                   fmt='.3f', desc=f'{engine} {method} log-log compression slope over spiked levels', between=(0, 2))
            st.add(f'tmt.{method}.{engine}.psms', row['quantified']['all']['psms'],
                   fmt=',', desc=f'{engine} {method} quantified PSMs at one percent', sign='+')
            st.add(f'tmt.{method}.{engine}.missing', 100 * row['missing_channels']['all_species_excl_126'],
                   fmt='.2f', desc=f'{engine} {method} percent missing reporter channels, 126 excluded', between=(0, 100))
