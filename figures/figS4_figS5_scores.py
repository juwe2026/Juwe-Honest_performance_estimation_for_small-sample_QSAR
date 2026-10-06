"""Fig. S4 and Fig. S5: what the Route A1 top-15 models do on the training data.

Fig. S4  one-component (LV1) apparent class separation, box plot with every compound overlaid
Fig. S5  two-component PLS-DA score plots with 95 % confidence ellipses
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from figstyle import *
import dataio as io

rcparams()
WB = io.xlsx("20260802_FigS5_S6_Rohdaten.xlsx")
RNG = np.random.default_rng(20260819)   # jitter only, fixed seed for reproducibility

# The source workbook stores the full strain name and the class as text.
SHEET_STRAIN = {"H.I.": "H. influenzae", "S.A.": "S. aureus", "S.Pneu": "S. pneumoniae",
                "S.Pyo": "S. pyogenes", "P.A.": "P. aeruginosa"}


def load(sheet):
    hdr, rows = io.sheet_rows(WB, sheet)
    return hdr, rows


def figS4(outdir):
    hdr, rows = load("FigS4_LV1_scores")
    # Sized to fit 180 mm width at 300 dpi; previously 276 mm wide.
    fig, axes = plt.subplots(1, 5, figsize=(7.15, 2.02), sharey=False)
    for c, s in enumerate(STRAINS):
        ax = axes[c]
        sub = [r for r in rows if r[0] == SHEET_STRAIN[s]]
        if not sub:
            raise SystemExit("no rows for %s in FigS4_LV1_scores" % s)
        groups = {}
        for r in sub:
            groups.setdefault(str(r[2]).strip(), []).append(float(r[3]))
        data = [groups.get("non-active", []), groups.get("active", [])]
        bp = ax.boxplot(data, positions=[0, 1], widths=0.55, showfliers=False,
                        patch_artist=True, zorder=3)
        for i, box in enumerate(bp["boxes"]):
            box.set(facecolor=tint(CLASS["non-active" if i == 0 else "active"], 0.68),
                    edgecolor=CLASS["non-active" if i == 0 else "active"], linewidth=0.9)
        for part in ("whiskers", "caps", "medians"):
            for j, ln in enumerate(bp[part]):
                ln.set(color=CLASS["non-active"] if j < len(bp[part]) / 2 else CLASS["active"],
                       linewidth=0.9)
        for i, vals in enumerate(data):
            x = i + RNG.uniform(-0.16, 0.16, len(vals))
            ax.plot(x, vals, "o", markersize=3.4,
                    color=CLASS["non-active" if i == 0 else "active"],
                    markeredgecolor="white", markeredgewidth=0.4, alpha=0.9, zorder=5)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["non-\nactive", "active"], fontsize=8.5)
        ax.set_title(STRAIN_LABEL[s], style="italic", fontsize=8.5)
        tidy(ax)
        ax.tick_params(axis="y", labelsize=8.0)       # revision (October 2026, EW): 10.5 -> 8 pt
        ax.yaxis.set_major_locator(plt.MultipleLocator(5))   # same tick spacing as submitted
        if c == 0:
            ax.set_ylabel("LV1 score (apparent)", fontsize=8.5)   # 11.5 -> 8.5 pt
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.text(0.5, 0.005, "Route A1 top-15 descriptors; every compound is shown, "
             "the LV1 axis is oriented so that active compounds score higher",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "FigS4_LV1_class_separation", outdir)


def ellipse(ax, x, y, colour):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)):
        return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    chi2_95 = 5.991                       # two degrees of freedom
    w, h = 2 * np.sqrt(np.maximum(vals, 0) * chi2_95)
    ax.add_patch(Ellipse((np.mean(x), np.mean(y)), w, h, angle=theta,
                         facecolor=colour, alpha=0.12, edgecolor=colour,
                         linewidth=1.0, linestyle=(0, (4, 3)), zorder=2))


def figS5(outdir):
    hdr, rows = load("FigS5_scores_2comp")
    # Sized to fit 180 mm width at 300 dpi; previously 289 mm wide. Trimmed slightly
    # further (7.03 -> 6.97 in) because the tight bbox padding pushed the rendered
    # PDF to 180.6 mm, 0.6 mm over the print area (round 4 report, 1.3).
    fig, axes = plt.subplots(1, 5, figsize=(6.97, 1.77))
    for c, s in enumerate(STRAINS):
        ax = axes[c]
        sub = [r for r in rows if r[0] == SHEET_STRAIN[s]]
        if not sub:
            raise SystemExit("no rows for %s in FigS5_scores_2comp" % s)
        # Component1 in this workbook is already oriented so "active compounds
        # score higher" holds for all five strains, matching Fig. S4's LV1_score
        # convention (Additional file 12 corrected H. influenzae's sign at the data
        # level; flipping it again here would silently undo that fix).
        for name in ("non-active", "active"):
            xs = np.array([float(r[3]) for r in sub if str(r[2]).strip() == name])
            ys = np.array([float(r[4]) for r in sub if str(r[2]).strip() == name])
            ellipse(ax, xs, ys, CLASS[name])
            ax.plot(xs, ys, "o", markersize=4.2, color=CLASS[name], label=name,
                    markeredgecolor="white", markeredgewidth=0.5, zorder=4)
        ax.axhline(0, color=GRID, linewidth=0.7, zorder=1)
        ax.axvline(0, color=GRID, linewidth=0.7, zorder=1)
        ax.set_title(STRAIN_LABEL[s], style="italic", fontsize=8.5)
        # revision (October 2026, EW): axis labels 11.5 -> 8.5 pt, tick labels 10.5 -> 8 pt
        ax.set_xlabel("component 1", fontsize=8.5)
        ax.tick_params(labelsize=8.0)
        if c == 0:
            ax.set_ylabel("component 2", fontsize=8.5)
        tidy(ax, grid_axis=None)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0.04, 1, 0.88])
    fig.legend(handles[:2], labels[:2], loc="upper center", ncol=2, handlelength=1.2, fontsize=8.0)   # 10.5 -> 8 pt
    fig.text(0.5, 0.005, "Route A1 top-15 descriptors, apparent two-component model; "
             "dashed outlines are 95 % confidence ellipses",
             ha="center", fontsize=7.5, color=MUTED)
    save(fig, "FigS5_PLSDA_score_plots", outdir)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    figS4(out); figS5(out)
