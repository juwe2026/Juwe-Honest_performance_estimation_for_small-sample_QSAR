"""Fig. 6: honest predictive ability and descriptor-set size for all four routes."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
from figstyle import *
import dataio as io

rcparams()
fa = io.npy("final_all.npy")
rb12 = io.routeB12()

SPEC = [("A1", "A1_nq", "A1_p", "A1_n"),
        ("A1u", "A1u_nq", "A1u_p", "A1u_n"),
        ("A2", "A2_nq", "A2_p", "A2_n"),
        ("B", None, None, None)]


def main(outdir):
    # Gleiche Breitenverhaeltnisse, damit die Balken in beiden Feldern gleich breit erscheinen
    # Sized to fit 180 mm width at 300 dpi; previously 251 mm wide. Height and layout
    # (axis label size, panel proportions, a/b label size and offset) matched to
    # Fig. 7's panel a/b, which uses the same rotated-label, two-panel layout.
    # Revision (October 2026, EW): same canvas, margins and panel spacing as Fig. 7
    # (6.3 x 3.35 in; left 0.085, right 0.985, bottom 0.33, top 0.80, wspace 0.48), so that
    # both figures, placed at the same width, have their y-axes on the same vertical lines
    # and x-axes of equal length; font sizes as in Fig. 7. Layout only, data unchanged.
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.35), gridspec_kw={"width_ratios": [1, 1]})
    xs = np.arange(len(STRAINS))
    w = 0.2
    q_by_route, p_by_route = {}, {}
    for k, (route, qk, pk, nk) in enumerate(SPEC):
        q, p, n = [], [], []
        for s in STRAINS:
            if route == "B":
                q.append(rb12["nlv1"][s]["q2"]); p.append(rb12["nlv1"][s]["p_q2"]); n.append(12)
            else:
                q.append(float(fa[s][qk])); p.append(float(fa[s][pk])); n.append(int(fa[s][nk]))
        q_by_route[route], p_by_route[route] = q, p
        pos = xs + (k - 1.5) * w
        bars = axes[0].bar(pos, q, w * 0.86, label=ROUTE_LABEL[route],
                           color=ROUTE[route], edgecolor=ROUTE[route], linewidth=0.7, zorder=3)
        for i, b in enumerate(bars):
            if STRAINS[i] == ARTEFACT:
                b.set_hatch(ARTEFACT_HATCH); b.set_edgecolor(MUTED)
        bars2 = axes[1].bar(pos, n, w * 0.86, color=ROUTE[route],
                            edgecolor=ROUTE[route], linewidth=0.7, zorder=3)
        for i, b in enumerate(bars2):
            if STRAINS[i] == ARTEFACT:
                b.set_hatch(ARTEFACT_HATCH); b.set_edgecolor(MUTED)

    # Asterisks sit at one common height per strain (bar group), not at each bar's own
    # height (round 4 report, 4.4). The glyph size has come down twice: 13 -> 9 pt in round 5,
    # which turned touching stars into a 0.42-0.51 mm gap, and 9 -> 6.75 pt here. The group
    # stays the same width, because the positions follow the bars; only the glyphs narrow, so
    # the gap roughly doubles and four adjacent stars can be counted at a glance
    # (round 6 report, 3.1).
    routes = [spec[0] for spec in SPEC]
    for i, s in enumerate(STRAINS):
        if s == ARTEFACT:
            continue
        group_vals = [q_by_route[r][i] for r in routes]
        group_top = max(group_vals)
        star_y = group_top + (0.045 if group_top >= 0 else 0.02)
        for k, r in enumerate(routes):
            if p_by_route[r][i] < 0.05:
                pos_x = xs[i] + (k - 1.5) * w
                # 8 pt on the larger canvas: about the size the 6.75 pt glyph had in print
                axes[0].text(pos_x, star_y, "*", ha="center", va="bottom",
                             fontsize=8.0, color=INK, zorder=4)

    axes[0].axhline(0, color=MUTED, linewidth=0.7, zorder=2)
    axes[0].set_ylabel("nested $Q^2$", fontsize=9.5)
    axes[0].set_ylim(-0.35, 0.92)
    # Single line, one point smaller than panel (a)'s label, rather than the two-line
    # form that pushed up into the panel letter's position (round 5 report, 1.2).
    axes[1].set_ylabel("descriptors entering PLS-DA", fontsize=9.5)   # as in Fig. 7
    # Two-line wrapping (round 4) broke each rotated name into a diagonal pair that no
    # longer read as one label (round 5 report, 1.1): reverted to single-line names at
    # the same 30 degree rotation, which was already shallower than the pre-round-4
    # angle and needs no further change on its own.
    strain_labels = [STRAIN_LABEL[s] for s in STRAINS]
    for ax in axes:
        tidy(ax)
        ax.tick_params(axis="y", labelsize=8.0)   # revision (October 2026, EW): 10.5 -> 8 pt
        ax.set_xticks(xs)
        # Revision (October 2026, EW): strain names 10.5 -> 7.5 pt, about the size of the legend
        ax.set_xticklabels(strain_labels, style="italic", fontsize=8.0,   # as in Fig. 7
                           rotation=30, ha="right", va="top", rotation_mode="anchor")
    axes[0].yaxis.set_major_locator(plt.MultipleLocator(0.5))   # tick spacing as submitted
    axes[1].yaxis.set_major_locator(plt.MultipleLocator(20))
    for i, ax in enumerate(axes):
        ax.text(-0.32, 1.06, "ab"[i], transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.33, top=0.80, wspace=0.48)   # = Fig. 7
    fig.legend(handles, labels, loc="upper center", ncol=2, columnspacing=1.2, handlelength=1.5, fontsize=7.5)
    fig.text(0.5, 0.005,
             "* marks p < 0.05 in the permutation test of the nested $Q^2$; hatched bars mark\n"
             "the P. aeruginosa artefact; Route B is fixed at 12 descriptors",
             ha="center", fontsize=7.2, color=MUTED, linespacing=1.4)   # as the notes of Fig. 7
    save(fig, "Fig6_order_comparison", outdir, xspan=XSPAN_FIG6_FIG7)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
