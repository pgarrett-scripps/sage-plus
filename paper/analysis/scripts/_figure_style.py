"""Consistent, print-sized styling for chapter figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INK = "#273444"
MUTED = "#596775"
GRID = "#E3E8ED"
TEAL = "#247D83"
PURPLE = "#806195"
COLORS = {"upstream": "#30688E", "plus": "#D57632"}
MARKERS = {"upstream": "o", "plus": "s"}
LINES = {"upstream": "-", "plus": (0, (4, 2))}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5, "text.color": INK,
    "axes.labelcolor": INK, "axes.labelsize": 10.5, "axes.labelpad": 7,
    "axes.titlesize": 11.5, "axes.titleweight": "bold", "axes.titlepad": 13,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#A2ABB3", "axes.linewidth": .7,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "legend.fontsize": 10, "legend.frameon": False,
    "legend.handlelength": 2.3, "legend.columnspacing": 1.8,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "svg.fonttype": "none", "svg.hashsalt": "sage-plus-chapter",
    "figure.constrained_layout.h_pad": .08,
    "figure.constrained_layout.w_pad": .08,
    "figure.constrained_layout.wspace": .08,
})


def panel(ax, letter, title, grid="y"):
    ax.set_title(f"{letter}  {title}", loc="left")
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=.65)
    ax.set_axisbelow(True)


def engine_style(engine, markers=True):
    result = dict(color=COLORS[engine], linestyle=LINES[engine], linewidth=1.7)
    if markers:
        result.update(marker=MARKERS[engine], markersize=5.5,
                      markerfacecolor="white" if engine == "plus" else COLORS[engine],
                      markeredgewidth=1.2)
    return result


def engine_legend(fig, markers=True, extra=()):
    handles = [Line2D([], [], label=label, **engine_style(engine, markers))
               for engine, label in (("upstream", "Sage"), ("plus", "Sage Plus"))]
    handles.extend(extra)
    fig.legend(handles=handles, loc="outside upper center", ncol=len(handles))


def save_figure(fig, target):
    """Keep the chapter raster and a scalable companion from the same plot."""
    fig.savefig(target, dpi=360, bbox_inches="tight", pad_inches=.06,
                metadata={"Software": None})
    vector = target.parent / "vector" / target.with_suffix(".svg").name
    vector.parent.mkdir(exist_ok=True)
    fig.savefig(vector, bbox_inches="tight", pad_inches=.06,
                metadata={"Date": None, "Creator": None})
    plt.close(fig)
