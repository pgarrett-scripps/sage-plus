"""Read the retained mass offset evaluation summary."""
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
SUMMARY = '../benchmarks/scientific-results/mass-offset-20260919/summary.json'
MASS_OFFSET_INPUTS = [SUMMARY]
LABEL = {'indexed': 'Indexed', 'offsets-1': 'One offset',
         'offsets-2': 'Two offsets', 'offsets-3': 'Three offsets'}


def load():
    return json.loads((PAPER / SUMMARY).read_text())


def scale_rows(summary):
    return {row['configuration']: row for row in summary['scale']}


def add_mass_offset_stats(st):
    summary = load()
    rows = scale_rows(summary)
    indexed, single = rows['indexed'], rows['offsets-1']

    # Only what the prose reads. Every configuration's complete values are in
    # the supporting table, which reads the summary directly.
    for configuration, row in rows.items():
        key = configuration.replace('-', '.')
        st.add(f'offset.{key}.search', row['search_stage_seconds'], fmt='.1f',
               desc=f'{LABEL[configuration]} median search stage seconds', sign='+')
    for configuration in ('indexed', 'offsets-1', 'offsets-3'):
        st.add(f"offset.{configuration.replace('-', '.')}.rss", rows[configuration]['peak_rss_mib'],
               fmt=',.0f', desc=f'{LABEL[configuration]} median peak resident MiB', sign='+')
    for configuration in ('indexed', 'offsets-3'):
        st.add(f"offset.{configuration.replace('-', '.')}.psms", rows[configuration]['psms'],
               fmt=',', desc=f'{LABEL[configuration]} accepted PSMs at one percent', sign='+')

    st.add('offset.peptides.indexed', indexed['database_peptides'], fmt=',',
           desc='Indexed search database peptides', sign='+')
    st.add('offset.peptides.offset', single['database_peptides'], fmt=',',
           desc='Offset search database peptides', sign='+')
    st.add('offset.peptides.reduction',
           100 * (1 - single['database_peptides'] / indexed['database_peptides']),
           fmt='.1f', desc='Offset reduction in indexed peptides percent', sign='+', between=(0, 100))
    st.add('offset.rss.reduction',
           100 * (1 - single['peak_rss_mib'] / indexed['peak_rss_mib']),
           fmt='.1f', desc='Offset reduction in peak resident memory percent', sign='+', between=(0, 100))
    st.add('offset.search.multiplier',
           single['search_stage_seconds'] / indexed['search_stage_seconds'],
           fmt='.1f', unit='times', desc='Search stage multiplier for one offset', sign='+', between=(0, 20))
    st.add('offset.search.per_offset',
           (rows['offsets-3']['search_stage_seconds'] - indexed['search_stage_seconds']) / 3,
           fmt='.1f', desc='Mean added search seconds per configured offset', sign='+', between=(0, 60))

    total_shared = sum(row['agreement']['shared_spectra'] for row in summary['localization'])
    total_same = sum(row['agreement']['same_peptidoform'] for row in summary['localization'])
    st.add('offset.agreement.spectra', total_shared, fmt=',',
           desc='Phosphopeptide spectra compared between search modes', sign='+')
    st.add('offset.agreement.same', total_same, fmt=',',
           desc='Spectra whose accepted peptidoform is identical in both modes', sign='+')
    st.add('offset.agreement.differing', total_shared - total_same, fmt=',',
           desc='Spectra whose accepted peptidoform differs between modes', between=(0, 100))

    for mode in ('indexed', 'offset'):
        correct = sum(row[mode]['correct_site_events'] for row in summary['localization'])
        incorrect = sum(row[mode]['incorrect_site_events'] for row in summary['localization'])
        st.add(f'offset.sites.{mode}.correct', correct, fmt=',',
               desc=f'{mode} localized site events consistent with synthesis', sign='+')
        st.add(f'offset.sites.{mode}.error', 100 * incorrect / (correct + incorrect), fmt='.2f',
               desc=f'{mode} inconsistent site event percent', sign='+', between=(0, 100))
        threshold = next(r for r in summary['entrapment'][mode] if r['nominal_q'] == 0.01)
        st.add(f'offset.entrapment.{mode}.fdp', 100 * threshold['combined_fdp'], fmt='.2f',
               desc=f'{mode} combined entrapment FDP percent at one percent peptide q',
               sign='+', between=(0, 100))
        st.add(f'offset.entrapment.{mode}.peptides', threshold['targets'] + threshold['entrapments'],
               fmt=',', desc=f'{mode} accepted paired peptides at one percent', sign='+')
