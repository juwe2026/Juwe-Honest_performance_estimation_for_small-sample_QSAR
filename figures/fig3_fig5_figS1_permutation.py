"""Fig. 3, Fig. 5 and Fig. S1: fully nested permutation tests.

Fig. 3  Route A1, full descriptor sets, all five strains   (perm_alt_*)
Fig. 5  Route A1 restricted to the top-15 by effect size   (perm_a1e_*)
Fig. S1 the H. influenzae panel of Fig. 3, enlarged
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from figstyle import *
import dataio as io

rcparams()


def _fit_fontsize(s, max_width_in, max_fs=8.3, min_fs=6.0, bold=False):
    """Largest font size (down to min_fs) at which s renders no wider than
    max_width_in, so the '<n> of 2,000 permutations below axis' annotation
    always fits the width of the panel it sits above, regardless of column
    count."""
    prop = FontProperties(family="DejaVu Sans", weight="bold" if bold else "normal")
    fs = max_fs
    while fs > min_fs:
        w_in = TextPath((0, 0), s, size=fs, prop=prop).get_extents().width / 72.0
        if w_in <= max_width_in:
            return fs
        fs -= 0.2
    return min_fs


def panel(ax, null, obs, xlabel, route_colour, tick_fs=None, label_fs=None,
          note_max_fs=8.3, note_y=1.03, labelpad=None):
    null = np.asarray(null, dtype=float)
    finite = null[np.isfinite(null)]
    lo = np.percentile(finite, 0.5)
    hi = max(np.percentile(finite, 99.9), obs)
    span = hi - lo
    lo, hi = lo - 0.05 * span, hi + 0.12 * span
    clipped = np.clip(finite, lo, hi)
    n_below = int((finite < lo).sum())
    ax.hist(clipped, bins=45, range=(lo, hi), color=NULL_FILL,
            edgecolor=NULL_EDGE, linewidth=0.4, zorder=3)
    p95 = np.percentile(finite, 95)
    ax.axvline(p95, color=PCT95, linestyle=(0, (4, 3)), linewidth=1.2, zorder=4)
    ax.axvline(obs, color=route_colour, linewidth=2.0, zorder=5)
    tidy(ax, grid_axis="y")
    ax.set_xlabel(xlabel) if label_fs is None else ax.set_xlabel(xlabel, fontsize=label_fs)
    if tick_fs is not None:
        ax.tick_params(axis="both", labelsize=tick_fs)
    if labelpad is not None:
        ax.xaxis.labelpad = labelpad
        ax.yaxis.labelpad = labelpad
    ax.set_xlim(lo, hi)
    if n_below:
        msg = "%d of 2,000 permutations below axis" % n_below
        ax_width_in = ax.get_position().width * ax.figure.get_figwidth()
        fs = _fit_fontsize(msg, ax_width_in * 0.98, max_fs=note_max_fs)
        ax.text(0.015, note_y, msg, transform=ax.transAxes, fontsize=fs,
                color=MUTED, va="bottom")


def draw(prefix, strains, stem, outdir, route_colour, figsize,
         tick_fs=None, label_fs=None, note_max_fs=8.3, name_fs=9,
         name_gap_pt=3.0, h_pad=None, legend_fs=None, labelpad=None, rect_top=None):
    """One permutation figure.

    name_gap_pt is the vertical gap, in points, between the baseline of the grey
    "<n> of 2,000 permutations below axis" note and the bottom of the bold strain
    name above it. Giving it in points rather than as a fraction of the axes keeps
    the gap constant when the panel height changes, which it does between the
    five-strain figures and the single-strain Fig. S1.
    """
    multi = len(strains) > 1
    fig, axes = plt.subplots(len(strains), 2, figsize=figsize, squeeze=False)
    note_y = 1.03
    for r, s in enumerate(strains):
        d = io.perm(prefix, s)
        panel(axes[r][0], d["q2s"], d["obs_q2"], "null $Q^2$", route_colour,
              tick_fs, label_fs, note_max_fs, note_y, labelpad)
        panel(axes[r][1], d["aucs"], d["obs_auc"], "null AUC", route_colour,
              tick_fs, label_fs, note_max_fs, note_y, labelpad)
        if label_fs is None:
            axes[r][0].set_ylabel("permutations")
        else:
            axes[r][0].set_ylabel("permutations", fontsize=label_fs)
    if h_pad is None:
        h_pad = 2.0 if multi else 1.0
    # rect_top is the share of the figure height left free above the axes for the legend.
    # It also sets how far the first bold strain name sits below the legend, because the
    # name is anchored to the first axes.
    if rect_top is None:
        rect_top = 0.965 if multi else 0.88
    fig.tight_layout(rect=[0, 0, 1, rect_top], h_pad=h_pad)

    # Placed after tight_layout, when the axes have their final height, so the gap
    # to the note below can be set in points and comes out the same in every panel
    # (round 5 report, 3.4; round 7 request).
    for r, s in enumerate(strains):
        ax = axes[r][0]
        ax_h_pt = ax.get_position().height * fig.get_figheight() * 72.0
        note_top_pt = note_y * ax_h_pt + note_max_fs      # top of the note, in points above the axes bottom
        y = (note_top_pt + name_gap_pt) / ax_h_pt
        ax.text(0.0, y, STRAIN_LABEL[s], transform=ax.transAxes,
                fontsize=name_fs, fontweight="bold", style="italic", va="bottom")

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=route_colour, lw=2.0, label="observed value"),
               Line2D([], [], color=PCT95, lw=1.2, ls=(0, (4, 3)), label="95th percentile of the null"),
               plt.Rectangle((0, 0), 1, 1, fc=NULL_FILL, ec=NULL_EDGE, lw=0.4,
                             label="null distribution, 2,000 permutations")]
    fig.legend(handles=handles, loc="upper center",
               ncol=2 if multi else 3, columnspacing=1.2, handlelength=1.3,
               fontsize=(7.5 if multi else 9.0) if legend_fs is None else legend_fs)
    save(fig, stem, outdir)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    # Fig. 3 and Fig. 5 are drawn at exactly the text width of the article template
    # (6.268 in = 159.2 mm, A4 less the 25.4 mm margins) instead of the former 127 mm. Two
    # things follow. The histograms roughly double in area, because at the old width the two
    # columns left each one only 44 x 17 mm while the axis numbering and labels at the shared
    # default sizes (10.5 and 11.5 pt) ate nearly 60 % of the figure height. And the figure is
    # placed in the document at its own size rather than scaled down, so the type sizes set
    # here are the sizes that print: 8 pt on the axis numbers and 9 pt on the axis labels,
    # a little below the effective sizes in Fig. 2 and comfortably above the 7 pt that print
    # production treats as the floor.
    # rect_top 0.940 rather than the 0.965 used elsewhere: at 0.965 only 1.7 mm separated
    # the two-line legend from the first bold strain name, barely more than the 1.4 mm
    # between the legend's own two lines, so the name read as a third legend entry rather
    # than as the heading of the panel below it. At 0.940 the gap is 5.1 mm and the reading
    # order is unambiguous; it costs about 1 mm of height per panel.
    MULTI = dict(tick_fs=8.0, label_fs=9.0, note_max_fs=7.0, name_fs=9.0,
                 name_gap_pt=3.0, h_pad=1.3, legend_fs=8.0, labelpad=2.0, rect_top=0.940)
    TEXT_WIDTH_IN = 6.268
    draw("perm_alt", STRAINS, "Fig3_nested_permutation_RouteA1_all5", out, ROUTE["A1"],
         (TEXT_WIDTH_IN, 8.72), **MULTI)
    draw("perm_a1e", STRAINS, "Fig5_nested_permutation_RouteA1_top15", out, ROUTE["A1"],
         (TEXT_WIDTH_IN, 8.72), **MULTI)
    # Fig. S1 is the enlarged single-strain panel and keeps the larger typography of the
    # main text figures it is compared against; only the strain-name gap is now set the
    # same way as above.
    draw("perm_alt", ["H.I."], "FigS1_nested_permutation_HI", out, ROUTE["A1"], (6.47, 3.33))
