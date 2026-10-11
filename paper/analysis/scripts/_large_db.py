"""Read the retained large-database evaluation summary."""
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
SUMMARY = '../benchmarks/scientific-results/large-db-20260925/summary.json'
LARGE_DB_INPUTS = [SUMMARY]
MULTIPLE_ORDER = [0, '1x', '3x', '10x', '30x', '100x', 'full']
OUTCOME = {'completed': 'completed', 'refused': 'refused by preflight',
           'guard': 'stopped by memory guard', 'ceiling': 'allocation failed at ceiling',
           'failed': 'failed'}


def load():
    return json.loads((PAPER / SUMMARY).read_text())


def ordered(rows, prefilter=None):
    """Rows in database-size order, optionally for one prefilter mode."""
    keep = {row['multiple']: row for row in rows
            if prefilter is None or row['prefilter'] == prefilter}
    return [keep[key] for key in MULTIPLE_ORDER if key in keep]


def by_name(rows):
    return {row['name']: row for row in rows}


def completed(row):
    return row.get('outcome') == 'completed'


def outcome(row):
    return OUTCOME[row['outcome']]


def multiple_label(multiple):
    if multiple == 0:
        return 'Human only'
    if multiple == 'full':
        return 'Human + full catalog'
    return f'Human + catalog {multiple}'


def database_label(series, row):
    if series in ('scaling', 'narrow'):
        return multiple_label(row['multiple'])
    if series == 'six_frame':
        return {'annotated': 'Annotated proteomes',
                'six-frame': 'Six-frame genomes'}[row['database']]
    return f"{row['sample']}, sample metagenome"


def pairs(summary):
    """Completed (prefilter, unfiltered) row pairs."""
    found = []
    for series in ('scaling', 'six_frame', 'metaproteome'):
        rows = by_name(summary[series])
        for name, row in rows.items():
            other = rows.get(name.replace('-prefilter', '-full'))
            if name.endswith('-prefilter') and other and completed(row) and completed(other):
                found.append((row, other))
    return found


def sweep(summary, min_matched, max_peaks=None):
    return next(row for row in summary['threshold_sweep']
                if row['min_matched'] == min_matched and row['max_peaks'] == max_peaks)


def add_large_db_stats(st):
    summary = load()
    databases = summary['databases']
    subsets = databases['subsets']
    scaling = by_name(summary['scaling'])
    narrow = by_name(summary['narrow'])
    six = by_name(summary['six_frame'])
    meta = by_name(summary['metaproteome'])

    st.add('large.catalog.proteins', databases['catalog_proteins_kept'], fmt=',',
           desc='Cleaned gut catalog proteins', sign='+')
    st.add('large.catalog.residues', subsets['1000']['residues'] / 1e9, fmt='.2f',
           desc='Cleaned gut catalog residues, billions', sign='+')
    st.add('large.catalog.multiple', subsets['1000']['residues'] / databases['human_residues'],
           fmt='.0f', desc='Full catalog residues as a multiple of human residues', sign='+')
    st.add('large.human.residues', databases['human_residues'] / 1e6, fmt='.1f',
           desc='Human reference residues, millions', sign='+')

    human = scaling['hek-human-prefilter']
    st.add('large.human.psms', human['accepted_psms'], fmt=',',
           desc='Human-only accepted PSMs', sign='+')

    for multiple in ('1x', '3x', '10x', '30x'):
        row = scaling[f'hek-igc-{multiple}-prefilter']
        st.add(f'large.{multiple}.fdp', row['combined_fdp'], fmt='.2f',
               desc=f'{multiple} combined entrapment FDP percent', between=(0, 5))
    for multiple in ('10x',):
        row = scaling[f'hek-igc-{multiple}-prefilter']
        st.add(f'large.{multiple}.retained', 100 * row['reference_peptides_retained']
               / row['reference_peptides'], fmt='.1f',
               desc=f'{multiple} percent of human-only peptides still accepted', between=(0, 100))
    for multiple in ('3x', '10x', '30x'):
        row = scaling[f'hek-igc-{multiple}-prefilter']
        counts = row['prefilter_counts']
        st.add(f'large.{multiple}.streamed', counts['streamed'] / 1e6, fmt='.1f',
               desc=f'{multiple} peptides streamed through the prefilter, millions', sign='+')
        st.add(f'large.{multiple}.kept', counts['retained'] / 1e6, fmt='.1f',
               desc=f'{multiple} peptides kept by the prefilter, millions', sign='+')
        st.add(f'large.{multiple}.retention', 100 * counts['retained'] / counts['streamed'],
               fmt='.0f', desc=f'{multiple} prefilter retention percent', between=(0, 100))
        st.add(f'large.{multiple}.rss', row['peak_rss_gib'], fmt='.1f',
               desc=f'{multiple} prefilter peak RSS GiB', between=(0, 16))
        st.add(f'large.{multiple}.minutes', row['wall_seconds'] / 60, fmt='.1f',
               desc=f'{multiple} prefilter wall minutes', sign='+')
    st.add('large.10x.psms', scaling['hek-igc-10x-prefilter']['accepted_psms'], fmt=',',
           desc='10x accepted PSMs', sign='+')
    full3 = scaling['hek-igc-3x-full']
    st.add('large.3x.full.rss', full3['peak_rss_gib'], fmt='.1f',
           desc='3x unfiltered peak RSS GiB', between=(0, 16))
    st.add('large.3x.full.minutes', full3['wall_seconds'] / 60, fmt='.1f',
           desc='3x unfiltered wall minutes', sign='+')
    st.add('large.10x.full.need', scaling['hek-igc-10x-full']['refused_estimate_gib'],
           fmt='.1f', desc='GiB the 10x unfiltered preflight said it needed', sign='+')
    st.add('large.30x.psms', scaling['hek-igc-30x-prefilter']['accepted_psms'], fmt=',',
           desc='30x accepted PSMs', sign='+')
    st.add('large.100x.need', scaling['hek-igc-100x-prefilter']['refused_estimate_gib'],
           fmt='.1f', desc='GiB the 100x preflight said the unmodified digest needed', sign='+')

    for multiple in ('10x', '30x'):
        row = narrow[f'hek-igc-{multiple}-mono-prefilter']
        counts = row['prefilter_counts']
        st.add(f'large.mono.{multiple}.retention', 100 * counts['retained'] / counts['streamed'],
               fmt='.0f', desc=f'{multiple} monoisotopic prefilter retention percent',
               between=(0, 100))
        st.add(f'large.mono.{multiple}.rss', row['peak_rss_gib'], fmt='.1f',
               desc=f'{multiple} monoisotopic prefilter peak RSS GiB', between=(0, 16))
        st.add(f'large.mono.{multiple}.psms', row['accepted_psms'], fmt=',',
               desc=f'{multiple} monoisotopic accepted PSMs', sign='+')

    compared = pairs(summary)
    st.add('large.pairs.compared', len(compared), fmt=',',
           desc='Completed prefilter and unfiltered pairs', sign='+')
    st.add('large.pairs.min.shared', min(100 * pre['peptides_shared_with_full']
                                         / full['accepted_peptides'] for pre, full in compared),
           fmt='.1f', desc='Lowest percent of unfiltered accepted peptides also accepted with '
           'the prefilter, across pairs', between=(90, 100))
    st.add('large.pairs.max.psm.loss', max(100 * (1 - pre['accepted_psms'] / full['accepted_psms'])
                                           for pre, full in compared),
           fmt='.1f', desc='Largest percent drop in accepted PSMs with the prefilter, across pairs',
           between=(0, 5))

    exact, default = sweep(summary, 1), sweep(summary, 3)
    for label, row in (('exact', exact), ('default', default)):
        counts = row['prefilter_counts']
        st.add(f'large.sweep.{label}.retention', 100 * counts['retained'] / counts['streamed'],
               fmt='.0f', desc=f'10x prefilter retention percent, {label} threshold',
               between=(0, 100))
        st.add(f'large.sweep.{label}.rss', row['peak_rss_gib'], fmt='.1f',
               desc=f'10x prefilter peak RSS GiB, {label} threshold', between=(0, 16))
        st.add(f'large.sweep.{label}.fdp', row['combined_fdp'], fmt='.2f',
               desc=f'10x combined entrapment FDP percent, {label} threshold', between=(0, 5))
    st.add('large.sweep.psm.loss', 100 * (1 - default['accepted_psms'] / exact['accepted_psms']),
           fmt='.1f', desc='Percent fewer 10x PSMs at the default threshold than at one match',
           between=(0, 5))

    annotated, frames = six['lfq-annotated-full'], six['lfq-six-frame-full']
    st.add('large.six.annotated.peptides', annotated['database_peptides'] / 1e6, fmt='.1f',
           desc='Annotated database peptides, millions', sign='+')
    st.add('large.six.frame.peptides', frames['database_peptides'] / 1e6, fmt='.1f',
           desc='Six-frame database peptides, millions', sign='+')
    for label, row in (('annotated', annotated), ('frame', frames)):
        st.add(f'large.six.{label}.microbial', row['microbial_peptides'], fmt=',',
               desc=f'{label} database accepted microbial peptides', sign='+')
    st.add('large.six.shared', frames['microbial_shared_with_other_database'], fmt=',',
           desc='Microbial peptides accepted from both databases', sign='+')
    st.add('large.six.unannotated', frames['microbial_unannotated'], fmt=',',
           desc='Six-frame microbial peptides absent from the annotated proteomes', sign='+')
    st.add('large.six.annotated.share', 100 * frames['microbial_annotated']
           / frames['microbial_peptides'], fmt='.1f',
           desc='Percent of six-frame microbial peptides found in annotated proteomes',
           between=(90, 100))

    metagenome = summary['campi_databases']['gut-db2mg-human.fasta']['GUT_DB2MG.faa']
    st.add('large.campi.proteins', metagenome['kept'], fmt=',',
           desc='CAMPI sample-specific metagenome proteins', sign='+')
    st.add('large.campi.residues', metagenome['residues'] / 1e6, fmt='.0f',
           desc='CAMPI sample-specific metagenome residues, millions', sign='+')
    st.add('large.campi.F06.need', meta['campi-F06-full']['refused_estimate_gib'], fmt='.1f',
           desc='GiB the CAMPI F06 unfiltered preflight needed', sign='+')
    for sample in ('F06', 'F05'):
        row = meta[f'campi-{sample}-prefilter']
        st.add(f'large.campi.{sample}.psms', row['accepted_psms'], fmt=',',
               desc=f'CAMPI {sample} accepted PSMs', sign='+')
        st.add(f'large.campi.{sample}.microbial', row['microbial_peptides'], fmt=',',
               desc=f'CAMPI {sample} accepted microbial peptides', sign='+')
        st.add(f'large.campi.{sample}.rss', row['peak_rss_gib'], fmt='.1f',
               desc=f'CAMPI {sample} prefilter peak RSS GiB', between=(0, 16))
