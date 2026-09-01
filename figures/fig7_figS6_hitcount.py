"""Fig. 7 and Fig. S6: how the reduction scheme changes the null distribution of the
number of descriptors passing the univariate effect-size filter."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
from figstyle import *
import dataio as io

rcparams()
H = io.npy("hits_three_variants.npy")

SCHEMES = [("A2", "A2", 2208), ("A1", "A1s", 182), ("A1u", "A1u", 182)]
ALPHA = 0.05


def binom_sd(n_tests):
    return np.sqrt(n_tests * ALPHA * (1 - ALPHA))


def fig7(outdir):
    # Sized to fit 180 mm width at 300 dpi; previously 222 mm wide.
    # Height reduced from 3.38 in: about 10 mm between the axis labels and the two
    # explanatory lines at the foot was empty. The panel proportions are unchanged,
    # only the dead space is removed (round 6 report, 4.3).
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.00),
                             gridspec_kw={"width_ratios": [1, 1], "wspace": 0.48})
    xs = np.arange(len(STRAINS))
    w = 0.24
    for k, (route, key, n_tests) in enumerate(SCHEMES):
        infl = [np.std(H[s]["null_" + key], ddof=1) / binom_sd(n_tests) for s in STRAINS]
        bars = axes[0].bar(xs + (k - 1) * w, infl, w * 0.86, label=ROUTE_LABEL[route],
                           color=ROUTE[route], edgecolor=ROUTE[route], linewidth=0.7, zorder=3)
        for i, b in enumerate(bars):
            if STRAINS[i] == ARTEFACT:
                b.set_hatch(ARTEFACT_HATCH); b.set_edgecolor(MUTED)
    axes[0].axhline(1, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
    axes[0].set_ylabel("null SD / binomial SD", fontsize=9.5)
    # The hatched P. aeruginosa / Route A2 bar reaches 10.02; without an explicit top
    # the autoscaled axis ended at 10 and clipped it against the frame (round 4 report, 4.5).
    axes[0].set_ylim(top=11)
    tidy(axes[0])

    obs = [H[s]["obs_A1u"] for s in STRAINS]
    nulls = [np.asarray(H[s]["null_A1u"], dtype=float) for s in STRAINS]
    bp = axes[1].boxplot(nulls, positions=xs, widths=0.5, whis=(2.5, 97.5), showfliers=False,
                         patch_artist=True, zorder=3)
    for i, box in enumerate(bp["boxes"]):
        box.set(facecolor=NULL_FILL, edgecolor=NULL_EDGE, linewidth=0.8)
        if STRAINS[i] == ARTEFACT:
            box.set_hatch(ARTEFACT_HATCH)
    for part in ("whiskers", "caps", "medians"):
        for ln in bp[part]:
            ln.set(color=NULL_EDGE, linewidth=0.8)
    for i, (x, o, nl) in enumerate(zip(xs, obs, nulls)):
        above = o > np.percentile(nl, 97.5)
        axes[1].plot([x], [o], marker="D", markersize=6.5,
                     color=ROUTE["A1u"] if not above else "#C1272D",
                     markeredgecolor="white", markeredgewidth=0.8, zorder=5)
    axes[1].axhline(182 * ALPHA, color=MUTED, linewidth=0.9, linestyle=(0, (1, 2)), zorder=2)
    axes[1].set_ylabel("descriptors passing the filter", fontsize=9.5)
    tidy(axes[1])

    for i, ax in enumerate(axes):
        ax.set_xticks(xs)
        ax.set_xticklabels([STRAIN_LABEL[s] for s in STRAINS], style="italic", rotation=30, ha="right")
        ax.text(-0.32, 1.06, "ab"[i], transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
    from matplotlib.lines import Line2D
    h1, l1 = axes[0].get_legend_handles_labels()
    h2 = [Line2D([], [], marker="D", ls="none", color=ROUTE["A1u"], markersize=6.5, label="observed count, within the null range"),
          Line2D([], [], marker="D", ls="none", color="#C1272D", markersize=6.5, label="observed count, above the null range")]
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.369, top=0.752, wspace=0.48)
    fig.legend(h1 + h2, l1 + [h.get_label() for h in h2], loc="upper center",
               bbox_to_anchor=(0.5, 1.0),
               ncol=3, columnspacing=1.4, handlelength=1.5, fontsize=7.5)
    fig.text(0.26, 0.011,
             "a: inflation of the null distribution of\nthe hit count relative to independent tests.\n"
             "Dashed line means independent tests (ratio = 1).",
             ha="center", fontsize=7.2, color=MUTED, linespacing=1.4)
    fig.text(0.76, 0.011,
             "b: null range of the Route A1u hit count,\n2.5th to 97.5th percentile of 1,000 permutations.\n"
             "Dotted line means 5 % of 182 tests, the naive expectation.",
             ha="center", fontsize=7.2, color=MUTED, linespacing=1.4)
    save(fig, "Fig7_hitcount_inflation", outdir)


def figS6(outdir):
    # Sized to fit 180 x 225 mm at 300 dpi. The width-only
    # scaling used in the first pass left too little height for a 3 x 5 grid with row
    # and column labels, which then overlapped; the 225 mm height budget was not
    # actually used, so the figure is taller here instead of just narrower.
    fig, axes = plt.subplots(3, 5, figsize=(7.0, 5.6))
    for r, (route, key, n_tests) in enumerate(SCHEMES):
        for c, s in enumerate(STRAINS):
            ax = axes[r][c]
            nl = np.asarray(H[s]["null_" + key], dtype=float)
            o = H[s]["obs_" + key]
            hi = max(nl.max(), o) * 1.06
            ax.hist(nl, bins=35, range=(0, hi), color=NULL_FILL, edgecolor=NULL_EDGE,
                    linewidth=0.3, zorder=3)
            ax.axvline(n_tests * ALPHA, color=MUTED, linestyle=(0, (1, 2)), linewidth=1.0, zorder=4)
            ax.axvline(np.percentile(nl, 95), color=PCT95, linestyle=(0, (4, 3)), linewidth=1.0, zorder=4)
            ax.axvline(o, color=ROUTE[route], linewidth=2.0, zorder=5)
            tidy(ax)
            ax.set_xlim(0, hi)
            ax.tick_params(labelsize=7.5)
            if r == 0:
                ax.set_title(STRAIN_LABEL[s], style="italic", fontsize=8.5)
            if c == 0:
                ax.set_ylabel("permutations", fontsize=8.5, labelpad=2)
                ax.text(-0.68, 1.22, "abc"[r], transform=ax.transAxes, fontsize=11,
                        fontweight="bold", va="top")
                ax.text(-0.85, 0.5, ROUTE_LABEL[route] + "\n(%d tests)" % n_tests,
                        transform=ax.transAxes, fontsize=8.0, fontweight="bold",
                        rotation=90, ha="center", va="center", color=ROUTE[route])
            if r == 2:
                ax.set_xlabel("descriptors passing", fontsize=8.0)
    from matplotlib.lines import Line2D
    handles = [plt.Rectangle((0, 0), 1, 1, fc=NULL_FILL, ec=NULL_EDGE, lw=0.3,
                             label="null distribution, 1,000 permutations"),
               Line2D([], [], color=MUTED, ls=(0, (1, 2)), lw=1.0, label="naive expectation, 5 % of the tests"),
               Line2D([], [], color=PCT95, ls=(0, (4, 3)), lw=1.0, label="95th percentile of the null"),
               Line2D([], [], color=INK, lw=2.0, label="observed count (coloured by route)")]
    fig.tight_layout(rect=[0.075, 0, 1, 0.92], h_pad=1.6, w_pad=0.6)
    fig.legend(handles=handles, loc="upper center", ncol=2, columnspacing=1.0,
               handlelength=1.2, fontsize=7.5)
    save(fig, "FigS6_hitcount_null_distributions", outdir)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    fig7(out); figS6(out)
