"""Fig. 2 and Fig. 4: apparent against nested performance for the full descriptor
sets (Fig. 2) and for the n:p-constrained top-15 sets (Fig. 4)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
from figstyle import *
import dataio as io

rcparams()
fa = io.npy("final_all.npy")
a2t = io.npy("a2_top15_results.npy")
rb12 = io.routeB12()          # 12 descriptors, C_Count removed (perfect collinearity)


def series(top15):
    """Six bars per strain: A1 apparent/nested, A2 apparent/nested, B 1 LV / 3 LV."""
    out = []
    for s in STRAINS:
        if top15:
            a1 = (fa[s]["A1e15_appq"], fa[s]["A1e15_nq"], fa[s]["A1e15_p"],
                  fa[s]["A1e15_appa"], fa[s]["A1e15_na"], fa[s]["A1e15_pa"])
            a2 = (a2t[s]["appq"], a2t[s]["nq"], a2t[s]["pq"],
                  a2t[s]["appa"], a2t[s]["na"], a2t[s]["pa"])
        else:
            a1 = (fa[s]["A1_appq"], fa[s]["A1_nq"], fa[s]["A1_p"],
                  fa[s]["A1_appa"], fa[s]["A1_na"], fa[s]["A1_pa"])
            a2 = (fa[s]["A2_appq"], fa[s]["A2_nq"], fa[s]["A2_p"],
                  fa[s]["A2_appa"], fa[s]["A2_na"], fa[s]["A2_pa"])
        n1, n3 = rb12["nlv1"][s], rb12["nlv3"][s]
        b1 = (n1["q2"], n1["p_q2"], n1["auc"], n1["p_auc"])
        b3 = (n3["q2"], n3["p_q2"], n3["auc"], n3["p_auc"])
        out.append(dict(a1=a1, a2=a2, b1=b1, b3=b3))
    return out


B_DARK = "#8E4A72"          # same hue as Route B, one step darker, for the 3-LV model

BARS = [
    ("Route A1, apparent", ROUTE["A1"], True),
    ("Route A1, nested",   ROUTE["A1"], False),
    ("Route A2, apparent", ROUTE["A2"], True),
    ("Route A2, nested",   ROUTE["A2"], False),
    ("Route B, 1 LV",      ROUTE["B"],  False),
    ("Route B, 3 LV",      B_DARK,      False),
]


def draw(top15, stem, outdir):
    data = series(top15)
    # Sized to fit the 180 x 225 mm print area at 300 dpi;
    # previously 191 x 183 mm, just over width.
    fig, axes = plt.subplots(2, 1, figsize=(6.41, 5.65), sharex=True)
    w = 0.135
    xs = np.arange(len(STRAINS))
    for row, metric in enumerate(["Q2", "AUC"]):
        ax = axes[row]
        for k, (label, colour, style) in enumerate(BARS):
            vals, ps = [], []
            for d in data:
                if k in (0, 1):
                    appq, nq, pq, appa, na, pa = d["a1"]
                elif k in (2, 3):
                    appq, nq, pq, appa, na, pa = d["a2"]
                if k == 0:
                    vals.append(appq if metric == "Q2" else appa); ps.append(None)
                elif k == 1:
                    vals.append(nq if metric == "Q2" else na); ps.append(pq if metric == "Q2" else pa)
                elif k == 2:
                    vals.append(appq if metric == "Q2" else appa); ps.append(None)
                elif k == 3:
                    vals.append(nq if metric == "Q2" else na); ps.append(pq if metric == "Q2" else pa)
                elif k == 4:
                    q, pq_, a, pa_ = d["b1"]
                    vals.append(q if metric == "Q2" else a); ps.append(pq_ if metric == "Q2" else pa_)
                else:
                    q, pq_, a, pa_ = d["b3"]
                    vals.append(q if metric == "Q2" else a); ps.append(pq_ if metric == "Q2" else pa_)
            face = tint(colour) if style is True else colour
            pos = xs + (k - 2.5) * w
            bars = ax.bar(pos, vals, w * 0.88, label=label if row == 0 else None,
                          color=face, edgecolor=colour, linewidth=0.7, zorder=3)
            for i, (b, v, p) in enumerate(zip(bars, vals, ps)):
                if STRAINS[i] == ARTEFACT:
                    b.set_hatch(ARTEFACT_HATCH)
                    b.set_edgecolor(MUTED)
                elif p is not None and p < 0.05:
                    ax.text(b.get_x() + b.get_width() / 2, v + (0.015 if v >= 0 else -0.055),
                            "*", ha="center", va="bottom", fontsize=10, color=INK, zorder=4)
        ax.axhline(0, color=MUTED, linewidth=0.7, zorder=2)
        tidy(ax)
        ax.set_ylabel("cross-validated $Q^2$" if metric == "Q2" else "AUC")
        ax.set_xticks(xs)
        ax.set_xticklabels([STRAIN_LABEL[s] for s in STRAINS], style="italic", fontsize=8.5, rotation=12, ha="center")
        if metric == "Q2":
            ax.set_ylim(-0.3, 1.02)
        else:
            ax.axhline(0.5, color=MUTED, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
            ax.set_ylim(0.45, 1.06)
            ax.set_xlim(xs[0] - 0.5, xs[-1] + 0.85)
            ax.text(xs[-1] + 0.58, 0.5 + 0.018, "chance", fontsize=7, color=MUTED, va="bottom", ha="left")
        ax.text(-0.12, 1.0, "ab"[row], transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.055),
               ncol=2, columnspacing=1.2, handlelength=1.3, fontsize=8.0)
    fig.text(0.5, -0.045,
             "pale fill = apparent, solid fill = nested; darker purple is the Route B\n"
             "three-component model; * marks p < 0.05; hatched bars mark P. aeruginosa",
             ha="center", fontsize=7.2, color=MUTED, linespacing=1.5)
    fig.tight_layout()
    save(fig, stem, outdir)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    draw(False, "Fig2_three_way_comparison", out)
    draw(True, "Fig4_three_way_comparison_top15", out)
