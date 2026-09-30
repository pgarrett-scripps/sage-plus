"""Read the known-site localization evaluation."""
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
SUMMARY = 'analysis/data/localization.json'
LOCALIZATION_INPUTS = [SUMMARY]
PRIMARY = 0.01


def load():
    return json.loads((PAPER / SUMMARY).read_text())


def point(curve, cutoff=PRIMARY):
    return next(row for row in curve['points'] if row['cutoff'] == cutoff)


def modification(per_type, *names):
    """The per-type row whose modification label names one of `names`."""
    return next(row for label, row in per_type.items() if any(name in label for name in names))


def add_localization_stats(st):
    summary = load()
    curves = summary['known_sites']
    st.add('loc.evaluable', curves['1']['evaluable'], fmt=',',
           desc='Evaluable single-site known-library phospho localizations', sign='+')
    for margin in summary['margins']:
        row = point(curves[str(margin)])
        st.add(f'loc.m{margin}.passing', row['passing'], fmt=',',
               desc=f'Known-site localizations at 1% localization q, margin {margin}', sign='+')
        st.add(f'loc.m{margin}.false', row['false'], fmt=',',
               desc=f'Wrong known-site localizations at 1%, margin {margin}', between=(0, 100000))
        st.add(f'loc.m{margin}.true', 100 * row['true_flr'], fmt='.2f',
               desc=f'True known-site FLR percent at 1%, margin {margin}', between=(0, 100))
        st.add(f'loc.m{margin}.estimated', 100 * row['estimated_flr'], fmt='.2f',
               desc=f'Estimated FLR percent (largest passing q) at 1%, margin {margin}',
               between=(0, 1))
    strict, default = point(curves['2']), point(curves['1'])
    st.add('loc.m2.loss', 100 * (1 - strict['passing'] / default['passing']), fmt='.0f',
           desc='Percent fewer known-site localizations at 1% with margin 2 than margin 1',
           between=(0, 100))
    loose = point(curves['1'], 0.05)
    st.add('loc.m1.passing05', loose['passing'], fmt=',',
           desc='Known-site localizations at 5% localization q, margin 1', sign='+')
    st.add('loc.m1.true05', 100 * loose['true_flr'], fmt='.2f',
           desc='True known-site FLR percent at 5%, margin 1', between=(0, 100))
    phospho = modification(summary['pxd007058'], 'Phospho', '79.96')
    oxidation = modification(summary['pxd007058'], 'Oxidation', '15.99')
    for name, row in (('phospho', phospho), ('oxidation', oxidation)):
        st.add(f'loc.bio.{name}.q01', row['localizations_q01'], fmt=',',
               desc=f'PXD007058 {name} localizations at 1%', sign='+')
        st.add(f'loc.bio.{name}.certain', row['certain'], fmt=',',
               desc=f'PXD007058 {name} localizations reported certain', between=(0, 1e6))
        st.add(f'loc.bio.{name}.sites', row['protein_sites_q01'], fmt=',',
               desc=f'PXD007058 {name} protein sites at 1% localization and site q', sign='+')
    groups = summary['localized_peptide']
    st.add('loc.proforma.rows', groups['rows'], fmt=',',
           desc='PSM rows with a ProForma localized peptide', sign='+')
    st.add('loc.proforma.groups', 100 * groups['scored_groups'] / groups['rows'], fmt='.1f',
           desc='Percent of localized peptides written with a scored position group',
           between=(0, 100))
