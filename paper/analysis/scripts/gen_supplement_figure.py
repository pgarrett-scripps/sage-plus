"""Visual summaries of supplemental comparisons from the frozen evidence."""
from statistics import median
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import StrMethodFormatter
from _assets import record
from _figure_style import plt, COLORS, INK, MUTED, GRID, TEAL, PURPLE, engine_style, engine_legend, panel, save_figure
from _scientific import PAPER, INPUTS, ENGINE, load
from _report import report
from _matched_fdp import load as matched_data, MATCHED_INPUTS


def finish(fig, name, inputs, desc):
    target = PAPER / 'figures' / f'supplement-{name}.png'
    save_figure(fig, target)
    record(f'fig.supplement-{name}', str(target.relative_to(PAPER)), kind='figure', inputs=inputs, desc=desc)
    vector = target.parent / 'vector' / target.with_suffix('.svg').name
    record(f'fig.supplement-{name}-vector', str(vector.relative_to(PAPER)), kind='figure', inputs=inputs,
           desc=f'Scalable companion: {desc}')


def dots(ax, x, y, engine, **kwargs):
    ax.plot(x, y, **{**engine_style(engine), 'linestyle': 'none'}, clip_on=False, **kwargs)


def row_axis(ax, names, positions=None):
    if positions is None:
        positions = np.arange(len(names))
    ax.set_yticks(positions, names)
    ax.set_ylim(positions[-1] + .55, positions[0] - .55)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=6)


def point_legend(fig):
    fig.legend(handles=[Line2D([], [], label=ENGINE[e], **{**engine_style(e), 'linestyle': 'none'})
                        for e in ('upstream', 'plus')], loc='outside upper center', ncol=2)


def resources(local=False):
    pilot, _ = load()
    keys = ('standard', 'common-mods') if local else ('human', 'hye')
    names = ['Standard', 'Common modifications'] if local else ['HEK', 'Mixture']
    metrics = [('wall_seconds', 'Wall time', 'Seconds', '.2f'),
               ('peak_rss_mib', 'Peak memory', 'MiB', ',.1f')]
    if local:
        metrics.append(('target_psms', 'Accepted PSMs', 'Target PSMs', ',.0f'))
    fig, axes = plt.subplots(1, len(metrics), figsize=(7.5, 3.45), sharey=True, layout='constrained')
    for ax, letter, (metric, title, xlabel, fmt) in zip(axes, 'ABC', metrics):
        all_values = []
        for i, key in enumerate(keys):
            for engine, offset in (('upstream', -.16), ('plus', .16)):
                group = [r for r in pilot['jobs'] if r['engine'] == engine and
                         ((r['suite'] == 'local-paired' and r['workload'] == key and not r['warmup'])
                          if local else r['suite'].startswith(f'entrapment-{key}-'))]
                assert len(group) == (3 if local else 6) and all(r['status'] == 'complete' for r in group)
                values = [r[metric] for r in group]
                value = median(values)
                all_values.extend(values)
                y = i + offset
                ax.scatter(values, y + np.linspace(-.055, .055, len(values)), s=10, alpha=.3,
                           color=COLORS[engine], zorder=2)
                dots(ax, value, y, engine, zorder=3)
                ax.annotate(format(value, fmt), (max(values), y), xytext=(7, 0),
                            textcoords='offset points', va='center', fontsize=9, color=INK)
        ax.set_xlim(0, max(all_values) * (1.5 if local else 1.28))
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.locator_params(axis='x', nbins=3 if local else 4)
        ax.set_xlabel(xlabel)
        row_axis(ax, names)
        panel(ax, letter, title, grid='x')
    point_legend(fig)
    finish(fig, 'local' if local else 'entrapment-resources', INPUTS,
           'Contextual local timing, memory, and PSM yield' if local else 'Entrapment resource measurements across files and seeds')


def threshold_yield():
    pilot, _ = load()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), layout='constrained')
    for ax, study, title, letter in zip(axes, ('human', 'hye'), ('HEK', 'Mixture'), 'AB'):
        for engine in ('upstream', 'plus'):
            cells = [r for r in pilot['calibration'] if r['suite'].startswith(f'entrapment-{study}-')
                     and r['job'].endswith('-' + engine)]
            qs = [.001, .005, .01, .02, .05]
            values = [np.mean([p['targets'] for r in cells for p in r['thresholds'] if p['nominal_q'] == q]) for q in qs]
            ax.plot(np.array(qs) * 100, values, **engine_style(engine))
        ax.set_xscale('log')
        ax.set_xticks([.1, .5, 1, 2, 5], ['0.1', '0.5', '1', '2', '5'])
        ax.minorticks_off()
        ax.set_xlabel('Nominal peptide q-value (%)')
        ax.set_ylabel('Mean target peptides')
        ax.set_ylim(0, 15000 if study == 'human' else 50000)
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        panel(ax, letter, title)
    engine_legend(fig)
    finish(fig, 'threshold-yield', INPUTS, 'Mean target peptide yield across prespecified nominal q-value thresholds')


def matched_yield():
    data = matched_data()
    fig, axes = plt.subplots(1, 3, figsize=(7.7, 3.9), sharey=True, layout='constrained')
    names = ['HEK · Sage', 'HEK · Sage Plus', 'Mixture · Sage', 'Mixture · Sage Plus']
    positions = [0, .75, 2, 2.75]
    for i, (study, engine) in enumerate((s, e) for s in ('human', 'hye') for e in ('upstream', 'plus')):
        row = next(r for r in data['summaries'] if r['study'] == study)['engines'][engine]
        dots(axes[0], row['mean_targets'], positions[i], engine)
        axes[0].annotate(f"{row['mean_targets']:,.1f}", (row['mean_targets'], positions[i]),
                         xytext=(6, 0), textcoords='offset points', va='center', fontsize=9)
        for ax, metric, lo, hi in ((axes[1], 'paired_fdp', 'fdp_min', 'fdp_max'),
                                    (axes[2], 'nominal_q', 'q_min', 'q_max')):
            values = [100 * r['selected'][metric] for r in data['runs'] if r['study'] == study and r['engine'] == engine]
            ax.hlines(positions[i], 100 * row[lo], 100 * row[hi], color=COLORS[engine], linewidth=1.5)
            ax.vlines([100 * row[lo], 100 * row[hi]], positions[i] - .065, positions[i] + .065,
                      color=COLORS[engine], linewidth=1)
            ax.scatter(values, positions[i] + np.linspace(-.08, .08, len(values)), s=13,
                       color=COLORS[engine], zorder=3)
    axes[0].set_xlim(0, 58000)
    axes[0].set_xticks([0, 20000, 40000], ['0', '20,000', '40,000'])
    axes[1].set_xlim(.978, 1.003)
    axes[1].set_xticks([.98, .99, 1], ['0.98', '0.99', '1.00'])
    axes[1].axvline(1, color=MUTED, linestyle=':', linewidth=1)
    axes[2].set_xlim(.75, 1.22)
    axes[2].set_xticks([.8, 1, 1.2])
    for ax, title, xlabel, letter in zip(axes, ('Target yield', 'Achieved FDP', 'Selected threshold'),
                                        ('Mean target peptides', 'Paired FDP (%)', 'Peptide q-value (%)'), 'ABC'):
        row_axis(ax, names, positions)
        ax.set_xlabel(xlabel)
        panel(ax, letter, title, grid='x')
    finish(fig, 'matched-fdp', MATCHED_INPUTS, 'Exploratory peptide yield, achieved FDP, and selected thresholds at the one-percent ceiling')


def public_counts():
    rows = report('public')
    names = ['HEK 1', 'HEK 2', 'A Alpha', 'B Alpha', 'A Beta', 'B Beta']
    positions = np.array([0, 1, 2.5, 3.5, 4.5, 5.5])
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 4.5), sharey=True, layout='constrained')
    for ax, metric, title, letter in zip(axes, ('target_psms', 'target_peptidoforms', 'decoy_psms'),
                                        ('Target PSMs', 'Peptidoforms', 'Decoy PSMs'), 'ABC'):
        for engine, offset in (('upstream', -.13), ('plus', .13)):
            values = [r['engines'][engine]['0.01'][metric] for r in rows]
            dots(ax, values, positions + offset, engine)
        ax.set_xlim(left=0)
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.locator_params(axis='x', nbins=3)
        ax.set_xlabel('Accepted count')
        row_axis(ax, names, positions)
        panel(ax, letter, title, grid='x')
    point_legend(fig)
    finish(fig, 'public-counts', ['analysis/data/report-extension/public.json'], 'Absolute accepted PSM, peptidoform, and decoy PSM counts for public inputs')


def shared_lfq():
    rows = report('lfq')['shared_ratio_pairs']
    species = ('human', 'yeast', 'ecoli')
    names = [f"{label}\n(n={rows[sp]['pairs']:,})" for sp, label in zip(species, ('Human', 'Yeast', 'E. coli'))]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True, layout='constrained')
    for engine, offset in (('upstream', -.15), ('plus', .15)):
        values = [rows[sp]['median_absolute_error'][engine] for sp in species]
        dots(axes[0], values, np.arange(3) + offset, engine)
        for y, value in enumerate(values):
            axes[0].annotate(f'{value:.3f}', (value, y + offset), xytext=(7, 0),
                             textcoords='offset points', va='center', fontsize=9)
    values = [rows[sp]['median_absolute_ratio_difference'] for sp in species]
    axes[1].scatter(values, np.arange(3), color=TEAL, marker='D', s=30)
    for y, value in enumerate(values):
        axes[1].annotate(f'{value:.4f}', (value, y), xytext=(7, 0),
                         textcoords='offset points', va='center', fontsize=9)
    axes[0].set_xlim(0, .6)
    axes[0].set_xticks([0, .2, .4, .6])
    axes[1].set_xlim(0, .01)
    axes[1].set_xticks([0, .005, .01], ['0', '0.005', '0.010'])
    for ax, title, xlabel, letter in zip(axes, ('Error against expected ratio', 'Difference between engines'),
                                        ('Median absolute log₂ error', 'Median absolute log₂ ratio difference'), 'AB'):
        row_axis(ax, names)
        ax.set_xlabel(xlabel)
        panel(ax, letter, title, grid='x')
    point_legend(fig)
    finish(fig, 'lfq-shared', ['analysis/data/report-extension/lfq.json'], 'Ratio error and between-engine differences on identical accepted pairs')


def ptm_diagnostic():
    pilot, _ = load()
    rows = [r for r in pilot['ptm'] if r['suite'] == 'ptm' and r['job'].endswith('-all')]
    fig, ax = plt.subplots(figsize=(7.2, 2.8), layout='constrained')
    for i, row in enumerate(rows):
        good, bad = row['correct_site_events'], row['incorrect_site_events']
        ax.barh(i, good, height=.5, color=TEAL)
        ax.barh(i, bad, left=good, height=.5, color=COLORS['plus'])
        ax.text(good / 2, i, f'{good:,}', va='center', ha='center', color='white', fontsize=10)
        ax.annotate(f"{bad:,} ({100 * row['empirical_site_error_fraction']:.2f}%)",
                         (good + bad, i), xytext=(6, 0), textcoords='offset points', va='center', fontsize=9)
        joint = row['joint_psm_peptide_localization_1pct']
        assert joint['correct_site_events'] + joint['incorrect_site_events'] == 0
    ax.set_xlim(0, 1600)
    ax.set_xticks([0, 500, 1000, 1500])
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    row_axis(ax, [f"HCD {r['job'].split('-')[1]}" for r in rows])
    ax.set_xlabel('Synthesis-consistent and inconsistent site events')
    ax.grid(axis='x', color=GRID, linewidth=.65)
    ax.set_axisbelow(True)
    fig.legend(handles=[Patch(facecolor=TEAL, label='Synthesis-consistent'),
                        Patch(facecolor=COLORS['plus'], label='Synthesis-inconsistent')],
               loc='outside upper center', ncol=2)
    finish(fig, 'ptm-diagnostic', INPUTS, 'Secondary Sage Plus synthesis consistency before the primary peptide filter')


if __name__ == '__main__':
    resources()
    resources(local=True)
    threshold_yield()
    matched_yield()
    public_counts()
    shared_lfq()
    ptm_diagnostic()
