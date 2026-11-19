"""Fig. 6 of the revised manuscript: ChEMBL benchmark, results per series, with the data sets of the study marked
as required by the pre-specification (22 October 2026).
(a) optimism, apparent minus nested Q2, per series and route
(b) calibrated inflation factor of the univariate hit count per series and route (log scale)
(c) calibrated inflation factor of Route A2 against the redundancy ratio G/p (H3), G/p on a logarithmic axis
(d) post hoc: observed calibrated factor divided by the factor predicted from the rank correlations of the tested
    descriptors (finite-sample model, no fitted parameter, naeherung_endlich_n.py, variant V2), Routes A1u and A2,
    one row per data set, endpoints of the study data sets shown individually
In (c) each data set of this study is shown as the median over its endpoints (open marker). The range over the endpoints
(monoterpenoids 10.2-10.8, peptides 12.7-13.4) is smaller than the marker on this axis, so no range bar is drawn.
Style as figures/figstyle.py of the code archive (Okabe-Ito route colours, DejaVu Sans, 300 dpi).
Filled markers: the 16 natural ChEMBL series. Open markers: data sets of this study.
Inputs: Additional file 12 (results workbook of the main analysis, sheet Summary), Additional file 2 (HitCounts_all_data_sets,
Table1_new_all_arms, Main_Table_1) and Additional file 7 (Gaucher re-analysis, sheet Routes_A1_A1u_A2).
Optimism of the monoterpenoids and peptides: apparent minus nested Q2 from revision2/analysis/scheinbare_q2_rev1code.py
(with the scripts in code/ of this archive, which reproduce Table 1 and Additional file 2 exactly).
Changes for the BMC figure requirements (6a, approved by EW on 06.11.2026 at 16:24 Berlin): every text at least 7.5 pt,
reading aid with the final reference (Additional file 9, section B6). Data and drawing otherwise unchanged.
Change of 07.11.2026 (approved by EW, item F01 of the review of 06.11.2026): "Peptides" instead of "Octapeptides",
also in the data-set names of Additional file 2. Data and drawing otherwise unchanged.
Change of 07.11.2026 (approved by EW, variant F6-C): 172.8 mm instead of 193.1 mm high, so that figure and legend fit
into 225 mm. Panels a and b flatter, the gap above panels c and d 5 mm smaller, slightly more space between panels a
and b, the reading aid in one line ("Additional file 9, B6"). Data, font sizes and line widths unchanged.
Change of 18.11.2026 (approved by EW, items A4 and A5 of 17.11.2026): the key names the open markers "comparison data
sets (open)" instead of "data sets of this study (open)", as the legend of Fig. 6 does, and panel c has a linear y axis
(0 to 15) instead of a logarithmic one, which spreads the 16 series. The x axis of panel c stays logarithmic. Data and
drawing otherwise unchanged.
Usage: python fig6_benchmark.py <AF12.xlsx> <AF2.xlsx> <AF7.xlsx> <scheinbar.json> <endlich_n.json> <outstem>
  scheinbar.json  revision2/results/Scheinbare_Q2_Octapeptide_Monoterpenoide 05.11.2026 2059.json
  endlich_n.json  revision2/results/Inflation_Naeherung_endliches_n_posthoc 05.11.2026 2232.json
  Writes <outstem>.png (300 dpi) and <outstem>.pdf."""
import sys, re, numpy as np, openpyxl, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

ROUTE = {"A1": "#E69F00", "A1u": "#009E73", "A2": "#0072B2"}
LABEL = {"A1": "Route A1", "A1u": "Route A1u", "A2": "Route A2"}
MARK = {"A1": "o", "A1u": "s", "A2": "^"}
INK, MUTED, GRID = "#1A1A1A", "#5A5A5A", "#D8D8D8"
ORDER = ["A1", "A1u", "A2"]                      # left to right within each column
OFF = {"A1": -0.24, "A1u": 0.0, "A2": 0.24}
FMIN = 7.5                                       # smallest text (pt), BMC: legible at final size (6a, EW 06.11.2026)
R1 = {"DS03", "DS11", "DS15"}                    # series with numerically unstable predictions
plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "font.family": "DejaVu Sans", "font.size": 9,
    "axes.labelsize": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.edgecolor": MUTED,
    "axes.linewidth": 0.7, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8.5, "legend.frameon": False,
    "grid.color": GRID, "grid.linewidth": 0.6, "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42})

def sheet(path, name):
    return list(openpyxl.load_workbook(path, read_only=True, data_only=True)[name].iter_rows(values_only=True))

# ---------- benchmark
rows = sheet(sys.argv[1], "Summary")
h = [i for i, r in enumerate(rows) if r and r[0] == "Data set"][0]; H = list(rows[h])
S = [dict(zip(H, r)) for r in rows[h + 1:] if r and r[0] and str(r[0]).startswith("DS") and int(r[0][2:]) <= 16]
S.sort(key=lambda d: d["Data set"])
ds = [d["Data set"] for d in S]
opt = {"A1": [d["Optimism A1 (apparent - nested Q2)"] for d in S], "A1u": [d["Optimism A1u"] for d in S],
       "A2": [d["Optimism A2"] for d in S]}
fac = {r: [d[f"Calibrated inflation {r}"] for d in S] for r in ORDER}
gp = [d["G/p"] for d in S]
groups = [(d["Activity type / subset"]) for d in S]

# ---------- data sets of this study
hc = sheet(sys.argv[2], "HitCounts_all_data_sets"); hh = list(hc[0])
key = {"Route A2": "A2", "Route A1": "A1", "Route A1u": "A1u"}
study_fac = {"Monoterpenoids": {r: [] for r in ORDER}, "Peptides": {r: [] for r in ORDER}, "Gaucher": {r: [] for r in ORDER}}
for r in hc[1:]:
    if not r or r[0] is None: continue
    d = dict(zip(hh, r))
    if d["Data_set"].startswith("Monoterpenoids") and d["Endpoint"] != "P. aeruginosa":     # four evaluable endpoints
        study_fac["Monoterpenoids"][key[d["Route"]]].append(d["Inflation_vs_own_rate"])
    elif d["Data_set"].startswith("Peptides, main"):
        study_fac["Peptides"][key[d["Route"]]].append(d["Inflation_vs_own_rate"])
    elif d["Data_set"].startswith("Gaucher"):
        study_fac["Gaucher"][key[d["Route"]]].append(d["Inflation_vs_own_rate"])
study_gp = {"Monoterpenoids": 182 / 2208, "Peptides": 50 / 1275, "Gaucher": 428 / 590}
# optimism of the study data sets: apparent minus nested Q2, both recomputed with the scripts in code/
import json as _json
SB = _json.load(open(sys.argv[4]))
study_opt = {"Monoterpenoids": {rt: [SB["Monoterpenoide"][k][rt]["scheinbar"] - SB["Monoterpenoide"][k][rt]["genestet"]
                                      for k in ["HI", "SA", "SPneu", "SPyo"]] for rt in ORDER},
             "Peptides": {rt: [SB["Octapeptide"][k][rt]["optimismus"] for k in ["EC", "SA", "PA", "AB"]] for rt in ORDER}}
g13 = sheet(sys.argv[3], "Routes_A1_A1u_A2"); gh = list(g13[0])
study_opt["Gaucher"] = {r: [] for r in ORDER}
for r in g13[1:]:
    if r and r[0] in ORDER and r[1] == 1:
        d = dict(zip(gh, r)); study_opt["Gaucher"][r[0]].append(d["Apparent_Q2"] - d["Nested_Q2"])

# ---------- figure
# Height. FIG_H = 7.3, GAP_CUT = 5 and HEAD_CUT = 3.5 give the layout approved by EW on 07.11.2026 (variant F6-C).
# FIG_H = 8.1, GAP_CUT = 0, HEAD_CUT = 0, HSPACE = 0.22 and the two-line reading aid gave the layout of 06.11.2026.
# A smaller FIG_H shortens only panels a and b, GAP_CUT (mm) shortens the gap between panel b and panels c and d,
# HEAD_CUT (mm) the header (one line of reading aid instead of two). Panels c and d keep their absolute size.
import os
FIG_H = float(os.environ.get("FIG6_H_IN", "7.3"))
GAP_CUT = float(os.environ.get("FIG6_GAP_CUT_MM", "5")) / 25.4
HEAD_CUT = 3.5 / 25.4
HSPACE = 0.30
_cell = 8.1 * (0.872 - 0.11) / (2 + 0.40)          # approved: gridspec(2, 1, height_ratios=[2.05, 1.3], hspace=0.40, top=0.872)
HEAD, FOOT, BOT_H, GAP = (1 - 0.872) * 8.1 - HEAD_CUT, 0.11 * 8.1, 2 * _cell * 1.3 / 3.35, 0.40 * _cell - GAP_CUT
fig = plt.figure(figsize=(6.3, FIG_H))
top = fig.add_gridspec(2, 1, hspace=HSPACE, top=1 - HEAD / FIG_H, bottom=(FOOT + BOT_H + GAP) / FIG_H)
bot = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.36, top=(FOOT + BOT_H) / FIG_H, bottom=FOOT / FIG_H)
axa = fig.add_subplot(top[0]); axb = fig.add_subplot(top[1], sharex=axa)
axc = fig.add_subplot(bot[0]); axd = fig.add_subplot(bot[1])
x = np.arange(16)
XS = {"Monoterpenoids": 17.2, "Peptides": 18.4, "Gaucher": 19.6}

def tidy(ax, grid="y"):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if grid: ax.grid(True, axis=grid, zorder=0); ax.set_axisbelow(True)

from matplotlib.transforms import blended_transform_factory, ScaledTranslation
SEP = "#9A9A9A"
def frame(ax, labels):
    ax.axvline(16.25, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    for i in range(1, 4): ax.axvline(4 * i - 0.5, color=SEP, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(-0.7, 20.3)
    if not labels:
        plt.setp(ax.get_xticklabels(), visible=False); return
    labs = [f"{d[2:]}{'†' if d in R1 else ''}" for d in ds]
    ax.set_xticks(list(x) + list(XS.values()))
    ax.set_xticklabels(labs + ["Monoterp.", "Peptides", "Gaucher"], fontsize=7.5)
    for t in ax.get_xticklabels()[16:]: t.set_rotation(90); t.set_color(INK)
    base = blended_transform_factory(ax.transData, ax.transAxes)
    tr = base + ScaledTranslation(0, -17 / 72, ax.figure.dpi_scale_trans)       # bracket line below the tick labels
    for i, name in enumerate(["Enzyme IC50", "MIC", "Receptor Ki", "Cytotoxicity"]):
        x0, x1 = 4 * i - 0.35, 4 * i + 3.35
        ax.plot([x0, x0, x1, x1], [0.03, 0, 0, 0.03], transform=tr, color=MUTED, lw=0.8, clip_on=False,
                solid_capstyle="butt")
        ax.annotate(name, xy=(4 * i + 1.5, 0), xycoords=("data", "axes fraction"), xytext=(0, -20),
                    textcoords="offset points", ha="center", va="top", fontsize=7.5, color=INK)
    ax.annotate("ChEMBL series (DS)", xy=(7.5, 0), xycoords=("data", "axes fraction"), xytext=(0, -33),
                textcoords="offset points", ha="center", va="top", fontsize=8.5, color=INK)

# (a) optimism
lo, hi = -0.45, 1.0
for rt in ORDER:
    v = np.array(opt[rt], float); xx = x + OFF[rt]
    inside = (v >= lo) & (v <= hi)
    axa.plot(xx[inside], v[inside], ls="none", marker=MARK[rt], ms=4.6, color=ROUTE[rt], mec="white", mew=0.5, zorder=3)
    for xi, vi in zip(xx[~inside], v[~inside]):
        up = vi > hi; yy = hi if up else lo
        axa.plot([xi], [yy], ls="none", marker="^" if up else "v", ms=5, color=ROUTE[rt], mec="white", mew=0.5,
                 zorder=4, clip_on=False)
        left = OFF[rt] < 0
        txt = f"{vi:.1f}".replace("-", "−")
        if up:
            axa.annotate(txt, (xi, yy), xytext=(-3 if left else 3, 3), textcoords="offset points",
                         fontsize=FMIN, color=MUTED, va="bottom", ha="right" if left else "left")
        else:
            axa.annotate(txt, (xi, yy), xytext=(6, 1), textcoords="offset points",
                         fontsize=FMIN, color=MUTED, va="bottom", ha="left")
    for name in ["Monoterpenoids", "Peptides", "Gaucher"]:
        vv = study_opt[name][rt]
        axa.plot([XS[name] + OFF[rt]] * len(vv), vv, ls="none", marker=MARK[rt], ms=4.6, mfc="white", mec=ROUTE[rt], mew=1.0, zorder=3)
axa.axhline(0, color=MUTED, lw=0.8)
axa.set_ylim(lo, hi); axa.set_ylabel("apparent − nested Q²")
axa.set_yticks([t for t in (-0.25, 0.0, 0.25, 0.5, 0.75, 1.0) if lo <= t <= hi])   # the steps of the approved figure, fixed
axa.set_title("a  Optimism of selection outside the validation (H1)", loc="left", pad=15)
tidy(axa); frame(axa, False)

# (b) calibrated inflation factor
for rt in ORDER:
    axb.plot(x + OFF[rt], fac[rt], ls="none", marker=MARK[rt], ms=4.6, color=ROUTE[rt], mec="white", mew=0.5, zorder=3)
    for name in XS:
        vv = study_fac[name][rt]
        axb.plot([XS[name] + OFF[rt]] * len(vv), vv, ls="none", marker=MARK[rt], ms=4.6, mfc="white", mec=ROUTE[rt], mew=1.0, zorder=3)
axb.set_yscale("log"); axb.set_ylim(0.9, 18)
axb.set_yticks([1, 1.5, 2, 3, 5, 10, 15]); axb.set_yticklabels(["1", "1.5", "2", "3", "5", "10", "15"]); axb.minorticks_off()
axb.axhline(1, color=MUTED, lw=0.8, ls=(0, (4, 3)))
axb.set_ylabel("null SD / binomial SD")
axb.set_title("b  Inflation of the null hit count (H2)", loc="left")
tidy(axb); frame(axb, True)

# (c) H3: G/p against factor A2
axc.plot(gp, fac["A2"], ls="none", marker=MARK["A2"], ms=4.6, color=ROUTE["A2"], mec="white", mew=0.5, zorder=3)
def study_mark(ax, xv, vv, rt):
    vv = np.asarray(vv, float); med = float(np.median(vv))
    ax.plot([xv], [med], ls="none", marker=MARK[rt], ms=6.4, mfc="white", mec=ROUTE[rt], mew=1.3, zorder=6)
    return med
def leader(ax, name, xy, xytext, ha):
    ax.annotate(name, xy=xy, xytext=xytext, textcoords="data", fontsize=FMIN, color=INK, ha=ha, va="center",
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.6, shrinkA=2, shrinkB=4.2), zorder=7)
labc = {"Monoterpenoids": ((0.115, 6.2), "left"), "Peptides": ((0.0165, 6.2), "left"), "Gaucher": ((0.30, 4.6), "center")}
for name in XS:
    med = study_mark(axc, study_gp[name], study_fac[name]["A2"], "A2")
    leader(axc, name, (study_gp[name], med), *labc[name])
rho = spearmanr(gp, fac["A2"]).statistic
axc.text(0.03, 0.04, f"16 ChEMBL series:\nSpearman ρ = {rho:.2f}".replace("-", "−"), transform=axc.transAxes, fontsize=FMIN, color=MUTED, va="bottom")
axc.set_xscale("log")
axc.set_xlim(0.015, 1.2); axc.set_ylim(0, 15)
axc.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1]); axc.set_xticklabels(["0.02", "0.05", "0.1", "0.2", "0.5", "1"])
axc.set_yticks([0, 5, 10, 15]); axc.set_yticklabels(["0", "5", "10", "15"]); axc.minorticks_off()
axc.set_xlabel("G/p (groups / descriptors)"); axc.set_ylabel("null SD / binomial SD")
axc.set_title("c  Redundancy, Route A2 (H3)", loc="left")
tidy(axc, grid="both")

# (d) post hoc: observed / predicted calibrated factor (finite-sample model V2, no fitted parameter), one row per data set
EN = _json.load(open(sys.argv[5]))
SN = {"Monoterpenoids": "Monoterpenoide", "Peptides": "Octapeptide", "Gaucher": "Gaucher"}
ROWS = list(range(16)) + [17.2, 18.4, 19.6]          # row positions, top to bottom (y axis inverted)
YOFF = {"A1u": -0.17, "A2": 0.17}
axd.axvspan(0.9, 1.1, color="#F0F0F0", lw=0, zorder=1)
axd.axvline(1, color=MUTED, lw=0.8, zorder=2)
for rt in ["A1u", "A2"]:
    r_ = np.array(fac[rt]) / np.array([EN["benchmark"][d][rt]["V2"] for d in ds])
    axd.plot(r_, np.arange(16) + YOFF[rt], ls="none", marker=MARK[rt], ms=4.4, color=ROUTE[rt], mec="white", mew=0.5, zorder=3)
    for name, yy in zip(SN, ROWS[16:]):
        vv = np.array(study_fac[name][rt]) / EN["studien"][SN[name]][rt]["V2"]
        axd.plot(vv, [yy + YOFF[rt]] * len(vv), ls="none", marker=MARK[rt], ms=4.4, mfc="white", mec=ROUTE[rt], mew=1.0, zorder=4)
for i in range(1, 4): axd.axhline(4 * i - 0.5, color=SEP, lw=0.7, ls=(0, (3, 2)), zorder=1)
axd.axhline(16.25, color=MUTED, lw=0.6, ls=(0, (2, 2)))
axd.set_ylim(20.2, -0.7)
axd.set_yticks(ROWS); axd.set_yticklabels([f"{d[2:]}{'†' if d in R1 else ''}" for d in ds] + ["Monoterp.", "Peptides", "Gaucher"], fontsize=FMIN)
axd.tick_params(axis="y", length=0, pad=2)
axd.set_xlim(0.84, 1.16); axd.set_xticks([0.9, 1.0, 1.1]); axd.set_xticklabels(["0.9", "1.0", "1.1"])
axd.set_xlabel("observed / predicted factor")
axd.set_title("d  Inflation predicted (post hoc)", loc="left")
tidy(axd, grid="x")

# legend at the top
hdl = [Line2D([], [], ls="none", marker=MARK[r], ms=5.5, color=ROUTE[r], mec="white", label=LABEL[r]) for r in ORDER]
hdl += [Line2D([], [], ls="none", marker="o", ms=5.5, color=MUTED, label="ChEMBL series (filled)"),
        Line2D([], [], ls="none", marker="o", ms=5.5, mfc="white", mec=MUTED, mew=1.0, label="comparison data sets (open)")]
fig.legend(handles=hdl, loc="upper center", bbox_to_anchor=(0.5, 1 - 0.005 * 8.1 / FIG_H), ncol=5, handletextpad=0.25, columnspacing=1.0, fontsize=8)
fig.text(0.5, 1 - 0.038 * 8.1 / FIG_H, "† numerically unstable predictions (Additional file 9, B6).   c: open marker median over endpoints.   d: grey band ±10 %.",
         ha="center", va="top", fontsize=FMIN, color=MUTED, linespacing=1.35)
# check: no two texts overlap (tick labels, axis titles, panel titles, annotations, key, reading aid)
fig.canvas.draw(); _r = fig.canvas.get_renderer()
_T = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
_B = [t.get_window_extent(_r) for t in _T]
_ov = [(_T[i].get_text()[:20], _T[j].get_text()[:20]) for i in range(len(_T)) for j in range(i + 1, len(_T))
       if _B[i].overlaps(_B[j]) and _B[i].intersection(_B[i], _B[j]).width > 0.5 and _B[i].intersection(_B[i], _B[j]).height > 0.5]
# the approved figure has 15 overlapping text boxes, all between neighbouring row labels of panel d (box height
# includes ascender and descender, the glyphs themselves do not touch). Any other overlap is an error.
_rowlab = lambda x: re.fullmatch(r"\d\d†?", x) is not None
_bad = [o for o in _ov if not (_rowlab(o[0]) and _rowlab(o[1]))]
print("text overlaps", len(_ov), "of them other than neighbouring row labels of panel d:", _bad)
assert not _bad and len(_ov) <= 15, _ov
stem = sys.argv[6]
for ext in ("png", "pdf"):
    fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.05)
print("rho", rho, "study_fac", {k: {r: [round(v, 2) for v in vv] for r, vv in d.items()} for k, d in study_fac.items()})
print("study_opt", {k: {r: [round(v, 3) for v in vv] for r, vv in d.items()} for k, d in study_opt.items()})
off = {r: [(ds[i], round(v, 2)) for i, v in enumerate(opt[r]) if v < lo or v > hi] for r in ORDER}; print("off-scale", off)
