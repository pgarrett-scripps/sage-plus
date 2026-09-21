"""Illustrate attachment identity and the two modification search paths."""
from matplotlib.patches import FancyBboxPatch
from _assets import record
from _figure_style import plt, INK, MUTED, TEAL, PURPLE, save_figure
from _scientific import PAPER


def box(ax, x, y, width, height, text, color):
    ax.add_patch(FancyBboxPatch((x, y), width, height,
        boxstyle='round,pad=0.012,rounding_size=0.025',
        linewidth=1, edgecolor=color, facecolor='white'))
    ax.text(x + width / 2, y + height / 2, text,
        ha='center', va='center', fontsize=10, color=color)


def arrow(ax, start, end, color=INK):
    ax.annotate('', xy=end, xytext=start,
        arrowprops=dict(arrowstyle='->', color=color, lw=1.4))


def main():
    fig = plt.figure(figsize=(7.3, 5.8))
    top = fig.add_axes((.03, .53, .94, .43))
    bottom = fig.add_axes((.03, .04, .94, .41))
    for ax in (top, bottom):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    top.text(0, .96, 'A  Attachment sites retain distinct identities',
             weight='bold', fontsize=12)
    residues = 'KSTGGKAPR'
    positions = [.20 + i * .078 for i in range(len(residues))]
    top.plot([.10, .91], [.49, .49], color='#BAC5CC', lw=2, zorder=1)
    top.text(.07, .49, 'N', ha='center', va='center', color=TEAL, fontsize=13,
             bbox=dict(boxstyle='circle,pad=.3', facecolor='white', edgecolor=TEAL))
    for i, (x, residue) in enumerate(zip(positions, residues)):
        color = PURPLE if residue == 'K' else MUTED
        top.text(x, .49, residue, ha='center', va='center', fontsize=14,
                 color=color, bbox=dict(facecolor='white', edgecolor='none', pad=2))
        top.text(x, .35, str(i + 1), ha='center', fontsize=9, color=MUTED)
    top.text(.94, .49, 'C', ha='center', va='center', color=MUTED, fontsize=13)
    top.annotate('Terminal group\npeptide_n_term:K', xy=(.07, .57), xytext=(.13, .77),
        ha='center', fontsize=10, color=TEAL,
        arrowprops=dict(arrowstyle='->', color=TEAL, lw=1.2))
    top.annotate('First residue\nfirst_residue:K', xy=(.20, .54), xytext=(.47, .77),
        ha='center', fontsize=10, color=PURPLE,
        arrowprops=dict(arrowstyle='->', color=PURPLE, lw=1.2))
    top.annotate('Internal residue\ninternal_residue:K', xy=(positions[5], .55),
        xytext=(.81, .77), ha='center', fontsize=10, color=PURPLE,
        arrowprops=dict(arrowstyle='->', color=PURPLE, lw=1.2))
    top.text(.5, .12, 'One named definition shares its occurrence limit across eligible sites.',
        ha='center', fontsize=10)
    bottom.text(0, 1.02, 'B  Search strategy and evidence reuse', weight='bold', fontsize=12)
    box(bottom, .02, .60, .29, .23, 'Indexed modification\nEnumerate before searching', TEAL)
    box(bottom, .64, .60, .33, .23, 'Mass offset modification\nPlace during searching', PURPLE)
    box(bottom, .30, .22, .40, .22, 'Scored peptidoform\nShared confidence assessment', INK)
    arrow(bottom, (.24, .59), (.41, .45), TEAL)
    arrow(bottom, (.76, .59), (.59, .45), PURPLE)
    arrow(bottom, (.50, .22), (.50, .12))
    bottom.text(.5, .07, 'Localization and typed PTM library', ha='center', fontsize=11)
    bottom.text(.5, -.08, 'Indistinguishable attachments do not become reusable site evidence.',
        ha='center', fontsize=9.5, color=MUTED)
    target = PAPER / 'figures/modification-model.png'
    save_figure(fig, target)
    for kind, path in (('fig.modification-model', target),
                       ('fig.modification-model-vector', target.parent / 'vector/modification-model.svg')):
        record(kind, str(path.relative_to(PAPER)), kind='figure',
            inputs=['analysis/data/development-changes.json'],
            desc='Attachment identity, modification search paths, and typed library reuse')


if __name__ == '__main__':
    main()
