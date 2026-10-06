"""Fig. S2 and Fig. S3: VIP values and VIP curves for the three descriptor sets.

Fig. S2  VIP per descriptor, apparent against mean nested VIP, one row per strain
         (a) Route A1 top-15   (b) Route B, 12 descriptors   (c) Route B, 25 descriptors
Fig. S3  the same values as curves with all five strains overlaid
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
from figstyle import *
import dataio as io

rcparams()
WB = io.xlsx("20260710_VIP_EffectSize_Tabellen.xlsx")

# Route B modelling set after removal of C_Count, which is perfectly collinear
# with EsterLacton_flag (r = rho = 1.000); see Materials and methods.
_VIPB = io.vip_routeB()
REDUCED_12 = list(_VIPB["D12"]["H.I."].keys())

STRAIN_COLOUR = {"H.I.": "#E69F00", "S.A.": "#0072B2", "S.Pneu": "#009E73",
                 "S.Pyo": "#CC79A7", "P.A.": "#8C8C8C"}


def chemdes(strain):
    hdr, rows = io.sheet_rows(WB, "VIP_ChemDes_" + strain)
    return [(r[0], float(r[1]), float(r[2])) for r in rows]


def lit(strain, reduced=False):
    """reduced=True gives the 12-descriptor Route B modelling set, recomputed;
    reduced=False the complete a priori set of 25, unchanged."""
    if reduced:
        d = _VIPB["D12"][strain]
        return [(k, v["apparent"], v["nested"]) for k, v in d.items()]
    hdr, rows = io.sheet_rows(WB, "VIP_Lit_" + strain)
    return [(r[0], float(r[3]), float(r[4])) for r in rows if r[3] is not None]


def panel(ax, getter, s, show_title=True, numbered=False):
    """One panel: horizontal VIP bars for one strain.

    numbered=True: the 15- and 25-descriptor columns no longer fit readable
    descriptor names once the figure is sized to print width without a
    further journal shrink, so rows are numbered 1..N in the same
    VIP-ascending order shown, and the names are given once in the caption /
    cross-referenced to the additional file that already lists them.
    """
    data = sorted(getter(s), key=lambda t: t[1])
    names = [d[0].split("::")[-1] for d in data]
    app = [d[1] for d in data]
    nst = [d[2] for d in data]
    y = np.arange(len(data))
    ax.barh(y, app, height=0.72, color=tint(ROUTE["A1"], 0.55), edgecolor=ROUTE["A1"],
            linewidth=0.6, label="apparent VIP", zorder=3)
    ax.plot(nst, y, marker="o", ls="none", markersize=4.2, color=INK,
            markeredgecolor="white", markeredgewidth=0.5, zorder=5, label="mean nested VIP")
    ax.axvline(1.0, color=MUTED, linestyle=(0, (3, 3)), linewidth=0.9, zorder=2)
    if numbered:
        # Every row still has its own bar; only every 5th row is numbered (plus the
        # first and last) so the tick labels themselves cannot overlap at this height.
        n = len(names)
        shown = sorted(set([0] + list(range(4, n, 5)) + [n - 1]))
        ax.set_yticks([y[i] for i in shown])
        ax.set_yticklabels([str(i + 1) for i in shown], fontsize=7.0)
        ax.tick_params(axis="y", length=2.5)
        print("    %s, %d descriptors, numbered 1-%d: " % (s, len(names), len(names))
              + ", ".join("%d=%s" % (i + 1, n) for i, n in enumerate(names)))
    else:
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=7.0)
    ax.set_ylim(-0.8, len(data) - 0.2)
    # enge x-Grenze: ohne den Vorgabe-Rand von 5 % nutzen die Balken die Flaeche aus
    ax.set_xlim(0, max(max(app), max(nst), 1.05) * 1.02)
    tidy(ax, grid_axis="x")
    ax.set_xlabel("VIP")
    if show_title:
        ax.set_title(STRAIN_LABEL[s], style="italic", fontsize=11.5)


def figS2(outdir):
    """Fuenf Zeilen fuer die Staemme, drei Spalten fuer die Deskriptorsaetze.

    Die fruehere Anordnung hatte fuenf Felder nebeneinander; die Deskriptornamen
    beanspruchten dann den halben Feldbereich und die Balken blieben kurz. Mit
    drei Spalten statt fuenf ist jedes Feld rund dreimal breiter, und da die
    Namensspalte gleich breit bleibt, werden die Balken etwa dreimal laenger.
    """
    getters = [chemdes, lambda s: lit(s, True), lambda s: lit(s, False)]
    titles = ["a  A1 top-15",
              "b  Route B, 12",
              "c  Route B, 25"]
    # Sized to fit 180 x 225 mm at 300 dpi; previously 365 x 443 mm.
    fig, axes = plt.subplots(5, 3, figsize=(6.9, 8.5))
    for r, s in enumerate(STRAINS):
        for c in range(3):
            panel(axes[r][c], getters[c], s, show_title=False, numbered=(c != 1))
        axes[r][0].text(-0.46, 0.5, STRAIN_LABEL[s], transform=axes[r][0].transAxes,
                        rotation=90, style="italic", fontsize=8, fontweight="bold",
                        ha="center", va="center", clip_on=False)
    for c in range(3):
        axes[0][c].set_title(titles[c], fontsize=8.5, fontweight="bold", pad=10)
    handles, labels = axes[0][0].get_legend_handles_labels()
    from matplotlib.lines import Line2D
    handles.append(Line2D([], [], color=MUTED, ls=(0, (3, 3)), lw=0.9, label="VIP = 1"))
    labels.append("VIP = 1")
    fig.tight_layout(rect=[0.02, 0, 1, 0.975])
    fig.legend(handles, labels, loc="upper center", ncol=3, handlelength=1.0, columnspacing=1.0, fontsize=6.8)
    save(fig, "FigS2_VIP_per_descriptor", outdir)


def figS3(outdir):
    # a, b und c untereinander: jedes Feld bekommt die volle Breite, die
    # Deskriptornamen auf der x-Achse werden dadurch lesbar
    # Sized to fit 180 x 225 mm at 300 dpi; previously 291 x 380 mm.
    fig, axes = plt.subplots(3, 1, figsize=(6.3, 7.3))
    specs = [(lambda s: lit(s, True), "Route B, 12 descriptors"),
             (lambda s: lit(s, False), "Route B, 25 descriptors"),
             (chemdes, "Route A1, top-15 (union of the five strain-specific sets)")]
    for i, (getter, title) in enumerate(specs):
        ax = axes[i]
        per = {s: dict((d[0], d[1]) for d in getter(s)) for s in STRAINS}
        # Sort key includes the descriptor name as a deterministic tie-break (round 5
        # report, cf. 2.1): without it, exact ties in mean VIP were ordered by Python's
        # set() iteration, which depends on the interpreter's hash seed and so is not
        # reproducible between runs -- confirmed to change Fig. S3c's panel-c axis order
        # for at least one tied pair (RDFM9 / MoRSEM18) across repeated runs.
        allnames = sorted(set().union(*[set(v) for v in per.values()]),
                          key=lambda n: (-np.mean([per[s].get(n, np.nan) for s in STRAINS
                                                   if n in per[s]]), n))
        x = np.arange(len(allnames))
        for s in STRAINS:
            y = [per[s].get(n, np.nan) for n in allnames]
            ax.plot(x, y, marker="o", markersize=5.0, linewidth=2.0,
                    color=STRAIN_COLOUR[s], label=STRAIN_LABEL[s],
                    linestyle="--" if s == ARTEFACT else "-", zorder=3)
        ax.axhline(1.0, color=MUTED, linestyle=(0, (3, 3)), linewidth=0.8, zorder=2)
        # Leave a margin below the smallest plotted VIP. Without it the lowest markers sat
        # on the axis line and were clipped in half, most visibly in panel (a)
        # (round 6 report, 4.3).
        _vals = [v for d in per.values() for v in d.values()]
        _lo, _hi = min(_vals), max(max(_vals), 1.0)
        _pad = 0.10 * (_hi - _lo)
        ax.set_ylim(_lo - _pad, _hi + _pad)
        if i == 2:
            # Panel c: ~60 descriptor names at full size
            # were unreadable once shrunk to print width. The panel's point is the
            # pattern of gaps between strains, not the individual names, so the axis
            # is numbered 1..N in the same rank order (every 5th number, plus the
            # first and last, so the labels themselves cannot overlap); the names
            # behind each number are the per-strain top-15 sets already listed in
            # Additional file 9.
            n = len(allnames)
            shown = sorted(set([0] + list(range(4, n, 5)) + [n - 1]))
            ax.set_xticks([x[j] for j in shown])
            ax.set_xticklabels([str(j + 1) for j in shown], rotation=0, fontsize=7.5)
            print("  Fig. S3c: %d descriptors, numbered 1-%d in this rank order:" % (n, n))
            print("   " + ", ".join("%d=%s" % (j + 1, nm.split("::")[-1]) for j, nm in enumerate(allnames)))
            # Export the position-to-descriptor-name key so it can be published as a
            # sheet in Additional file 8 (round 5 report, 3.2), instead of only being
            # printed to stdout. Keep the "Toolkit::Descriptor" full name (matching the
            # Descriptor column of the existing VIP_ChemDes_* sheets), not the
            # stripped short label used for the console summary above: two distinct
            # descriptors here share the short name "FPSA-3" (PaDEL_3D and
            # BlueDesc_mmff94), which would make the key ambiguous otherwise.
            with open(os.path.join(outdir, "FigS3c_axis_key.csv"), "w") as fh:
                fh.write("Position,Descriptor\n")
                for j, nm in enumerate(allnames):
                    fh.write("%d,%s\n" % (j + 1, nm))
        else:
            ax.set_xticks(x)
            ax.set_xticklabels([n.split("::")[-1] for n in allnames], rotation=90, fontsize=7)
        tidy(ax)
        ax.set_ylabel("VIP" if i == 0 else "")
        ax.set_title(title, fontsize=10, fontweight="bold", pad=8)
        ax.text(-0.07, 1.05, "abc"[i], transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="bottom")
    handles, labels = axes[0].get_legend_handles_labels()
    # h_pad was 5.0, which left an empty band of about 20 mm between (a) and (b) and
    # again between (b) and (c): the rotated tick labels reserve their own space through
    # tight_layout, so the padding was pure slack on a figure already close to the
    # 225 mm height limit (round 6 report, 4.3).
    fig.tight_layout(rect=[0, 0, 1, 0.96], h_pad=1.2)
    fig.legend(handles, labels, loc="upper center", ncol=3, handlelength=1.3, fontsize=8.0)
    save(fig, "FigS3_VIP_curves", outdir)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    figS2(out); figS3(out)
