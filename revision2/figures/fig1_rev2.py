# -*- coding: utf-8 -*-
"""
Change of 17.11.2026 (approved by EW, item 6 of 17.11.2026): the box of the univariate filter reads
"the |effect size| > 0.3" instead of "the |r_rb| > 0.3", whose subscript was printed as a literal underscore. The
effect size is the rank-biserial correlation defined in the Methods. fig1_rev2c.py, which draws the figure of the article,
reads this wording and was changed in the same way. Nothing else changed.

Fig. 1, Revision 2 (draft of 06.11.2026): as fig1_v3.py, but without the marks E1 to E9 (the list E1 to E9 was
removed from the Methods), without the panel on the two-stage VIP procedure (now in the Supplement), with short panel
titles and a new panel d on replication and generalisation (octapeptides, Gaucher, ChEMBL benchmark).
Version 06.11.2026 03:10 (EW): panel c with two separate model boxes, PLS-DA (1 LV, leave-one-out) and random forest
(300 trees, out-of-bag), because each of the three pools was modelled with both models (D_rf.json: all2208, 182 centrality,
Morgan2 bits). Text in the boxes "Same pipeline" and "Control arms" centred vertically (EW 04:07 SAST). Result boxes
in panel d of equal height, centres at 1/4 and 3/4 of the height of the box "Same pipeline" (EW 04:12 SAST).

Original header of fig1_v3.py: Fig. 1, version 04.10.2026 (request EW): print size following the Springer guidelines for figures.

Springer artwork guidelines: font 8-12 pt at print size, figure width at most 174 mm
(two columns), height at most 234 mm. The previous version (fig1_v2.py, replaced by this file) was 180 mm wide, with 6-7 pt fonts.
This version draws directly in millimetres on a 174 mm wide canvas, every font is >= 8 pt
(text 8 pt, box titles 8.5 pt bold, panel letters 10 pt). Font Liberation Sans (metrically equal to Arial).
Content, colours, arrows and layout as fig1_v2.py, changed only:
  - line breaks in the boxes (more lines, because the font is larger),
  - E5/E6/E7 as corner marks in the input boxes of panel d (like E1-E4), no longer on the arrows,
  - legend row on two lines (3 + 3 entries),
  - the grouping of the survivors under Route A2 green with an orange frame ("input selected with labels").
  - Panel d: the Morgan2 arm with nested univariate filter as a separate orange branch with its own PLS-DA box
    (sensitivity analysis, not a control arm, audit round 2, 04.10.2026).
Subject matter: Route A2 filters all 2,208 descriptors first and only then groups the survivors. Therefore
the arrow to Route A2 starts from the box of the 2,208 descriptors, and its grouping lies inside the loop.
Checks at the end: every text lies inside its box, texts of one box do not overlap, boxes and
panel letters do not overlap, every font >= 8 pt, width = 174 mm, height <= 234 mm.

Usage in the code archive:  cd figures && python ../revision2/figures/fig1_rev2.py output
  writes output/Fig1_scheme_revised.png (600 dpi) and output/Fig1_scheme_revised.pdf.
Without an argument the script writes the time-stamped working files of the revision instead."""
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)
import sys, os, datetime, zoneinfo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = 'Liberation Sans'
plt.rcParams['pdf.fonttype'] = 42          # TrueType embedded
plt.rcParams['svg.fonttype'] = 'path'

GREEN = '#2E7D32'; GREENF = '#E7F3E8'; ORANGE = '#D96A00'; ORANGEF = '#FCECD9'
BLUE = '#1A5FB4'; BLUEF = '#E3EBF8'; GREY = '#5A5A5A'; GREYF = '#EFEFEF'; LINE = '#505050'
FB, FT, FTAG, FN, FP = 8.0, 8.5, 8.0, 8.0, 10.0      # text, box titles, E marks, notes, panel letters (pt)
MIN_PT = 8.0
LS = 1.25                                             # line spacing
WMM, HMM = 174.0, 234.0                               # figure size in mm
PAD = 1.8                                             # inner margin at the top (mm)
LT = 4.3                                              # height of the title line (mm)

fig = plt.figure(figsize=(WMM / 25.4, HMM / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, WMM); ax.set_ylim(0, HMM); ax.axis('off')
def Y(top): return HMM - top                          # distance from the top (mm) -> y coordinate

# columns (mm): gap after column 0 wider because of the elbow arrows
C0, W0 = 2.0, 35.0
C1, W1 = 48.0, 38.0
C2, W2 = 93.0, 38.0
C3, W3 = 138.0, 34.0
BOX = []; TXT = []; LETTERS = []

def box(x, top, w, h, txt, ec, fc, bold=False, title=None, tag=None, tc=None, center=False):
    """Box with upper edge 'top' (mm from the top). tc: title colour, if it differs from the frame."""
    y = Y(top) - h
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.5', linewidth=0.9,
                       edgecolor=ec, facecolor=fc, zorder=3)
    ax.add_patch(p); BOX.append((p, txt[:20]))
    own = []
    if title:
        off = 0.0
        if center:                                    # title + text as a block, vertically centred
            blk = LT + txt.count('\n') * FB * LS * 25.4 / 72 + FB * 25.4 / 72
            off = max(0.0, (h - 2 * PAD - blk) / 2)
        top = top + off
        own.append(ax.text(x + w / 2, Y(top) - PAD, title, ha='center', va='top', fontsize=FT, fontweight='bold',
                           color=tc or ec, zorder=4))
        own.append(ax.text(x + w / 2, Y(top) - PAD - LT, txt, ha='center', va='top', fontsize=FB, color='#1A1A1A',
                           zorder=4, linespacing=LS))
    else:
        own.append(ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=FB, color='#1A1A1A', zorder=4,
                           fontweight='bold' if bold else 'normal', linespacing=LS))
    if tag:
        own.append(ax.text(x + w - 1.0, y + h - 0.9, tag, ha='right', va='top', fontsize=FTAG, color='#666666',
                           fontweight='bold', zorder=5))
    TXT.append((own, p))
    return p

def ar(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=8, linewidth=0.85, color=LINE,
                                 zorder=2, shrinkA=0, shrinkB=0))

def line(xs, ys):
    ax.add_line(Line2D(xs, ys, linewidth=0.85, color=LINE, zorder=2, solid_capstyle='butt', solid_joinstyle='miter'))

def elbow(x1, y1, xm, y2, x2):
    """horizontally to the right up to xm, vertically to y2, then horizontally with arrow up to x2"""
    line([x1, xm], [y1, y1]); line([xm, xm], [y1, y2]); ar(xm, y2, x2, y2)

def panel(top, t, title=''):
    LETTERS.append(ax.text(0.2, Y(top), t, fontsize=FP, fontweight='bold', color='#1A1A1A', ha='left', va='top', zorder=6))
    if title:
        LETTERS.append(ax.text(4.2, Y(top) - 0.35, title, fontsize=FT, fontweight='bold', color='#1A1A1A', ha='left', va='top', zorder=6))

def note(xc, top, t, color):
    return ax.text(xc, Y(top), t, ha='center', va='top', fontsize=FN, color=color, style='italic', linespacing=LS)

def mid(top, h): return Y(top + h / 2)

# ============================================ (a) Route A
panel(0.5, 'a', 'Route A, worked example')
FX0, FX1 = C1 - 1.8, C2 + W2 + 1.8
FT_, FB_ = 11.5, 69.5                                  # frame top/bottom (mm from the top)
ax.add_patch(Rectangle((FX0, Y(FB_)), FX1 - FX0, FB_ - FT_, linewidth=1.0, edgecolor=ORANGE, facecolor='none',
                       linestyle=(0, (4, 2.6)), zorder=1))
NOTES = [ax.text((FX0 + FX1) / 2, Y(FT_) + 1.0, 'selection steps repeated inside every leave-one-out\nfold and under each of the 2,000 label permutations',
                 ha='center', va='bottom', fontsize=FN, color=ORANGE, style='italic', linespacing=LS)]
tD, hD = 13.0, 9.5;  tP, hP = 27.0, 13.5;  tG, hG = 45.0, 21.0
box(C0, tD, W0, hD, '31 monoterpenoids,\n5 binary endpoints', GREY, 'white', bold=True)
box(C0, tP, W0, hP, '2,208 cleaned\ndescriptors from\n5 libraries', GREEN, GREENF)
box(C0, tG, W0, hG, '|r| or |rho| > 0.7\n2,208 → 182 groups\nties → largest variance', GREEN, GREENF, title='Correlation grouping')
ar(C0 + W0 / 2, Y(tD + hD), C0 + W0 / 2, Y(tP)); ar(C0 + W0 / 2, Y(tP + hP), C0 + W0 / 2, Y(tG))
NOTES.append(note(C0 + W0 / 2, tG + hG + 0.8, 'label-free, computed once', GREEN))
tA2, hA2 = 13.0, 18.0;  tA1, hA1 = 34.5, 14.5;  tA1u, hA1u = 52.5, 14.5
box(C1, tA2, W1, hA2, 'univariate filter\napplied to all\n2,208 descriptors', ORANGE, ORANGEF, title='Route A2')
box(C1, tA1, W1, hA1, 'representative per group\nchosen by effect size', ORANGE, ORANGEF, title='Route A1')
box(C1, tA1u, W1, hA1u, 'representative per group\nchosen by centrality', GREEN, GREENF, title='Route A1u')
yA2, yA1, yA1u = mid(tA2, hA2), mid(tA1, hA1), mid(tA1u, hA1u)
assert Y(tG + hG) < yA1u < Y(tG), 'A1u arrow must start from the grouping box'
# arrows without crossing: 2,208 descriptors -> Route A2 (filter first), grouping -> Route A1 and A1u
xA2m, xGm = C0 + W0 + 3.0, C0 + W0 + 6.5
elbow(C0 + W0, mid(tP, hP), xA2m, yA2, C1)
ar(C0 + W0, yA1u, C1, yA1u)
line([xGm, xGm], [yA1u, yA1]); ar(xGm, yA1, C1, yA1)
ax.plot([xGm], [yA1u], marker='o', markersize=2.2, color=LINE, zorder=3)
# grouping of the survivors: label-free step (green), input from the label-dependent filter (orange frame)
tS, hS = 13.0, 18.0;  tF, hF = 38.0, 18.0
box(C2, tS, W2, hS, 'of the survivors,\ninput selected\nwith labels', ORANGE, GREENF, title='Correlation grouping', tc=GREEN)
box(C2, tF, W2, hF, 'Mann–Whitney p ≤ 0.05,\nthe |effect size| > 0.3\nthreshold never binds', ORANGE, ORANGEF, title='Univariate filter')
assert Y(tS + hS) < yA2 < Y(tS)
ar(C1 + W1, yA2, C2, yA2)
ar(C1 + W1, yA1, C2, Y(tF + 5.0))
ar(C1 + W1, yA1u, C2, Y(tF + hF - 4.0))
tM, hM = 18.0, 14.5
box(C3, tM, W3, hM, 'PLS-DA\n1 latent variable\nleave-one-out', BLUE, BLUEF, bold=True)
assert Y(tM + hM) < yA2 < Y(tM)
ar(C2 + W2, yA2, C3, yA2)                              # grouping of the survivors -> model (horizontal)
# symmetric to the centre of the filter box: start +/- 2 mm, end +/- d
yFc = mid(tF, hF); d = Y(tM + hM - 2.5) - yFc       # > 0: end at the top inside the model box, symmetric at the bottom
hN = 11.0; tN = (HMM - (yFc - d)) - hN / 2; assert d > 0
box(C3, tN, W3, hN, 'Null distribution\nof the hit count', ORANGE, ORANGEF)
ar(C2 + W2, yFc + 2.0, C3, yFc + d)
ar(C2 + W2, yFc - 2.0, C3, yFc - d)
A_END = max(FB_, tN + hN, tG + hG + 5.0)

# ============================================ (b) Route B
tb = A_END + 3.5; panel(tb, 'b', 'Route B'); tb += 5.0
hb0, hb1, hb3 = 13.5, 14.5, 10.0; cb = tb + hb1 / 2
box(C0, cb - hb0 / 2, W0, hb0, '25 literature\ndescriptors\nreduced to 12', GREEN, GREENF)
box(C1, cb - hb1 / 2, W1, hb1, 'no supervised selection,\nso nothing to nest', GREEN, GREENF, title='Route B')
box(C3, cb - hb3 / 2, W3, hb3, 'PLS-DA\n1 and 3 LV', BLUE, BLUEF, bold=True)
ar(C0 + W0, Y(cb), C1, Y(cb)); ar(C1 + W1, Y(cb), C3, Y(cb))
B_END = cb + hb1 / 2

C_END = B_END

# ============================================ (d) control arms
td = C_END + 3.5; panel(td, 'c', 'Control arms and a sensitivity branch'); td += 5.0
rows = [(10.5, 'all 2,208\ndescriptors', ''), (14.0, '182 label-free\nrepresentatives\n(panel a)', ''),
        (10.5, 'Morgan2\nfingerprint bits', '')]
gap = 2.5; t = td; mids = []
for h, txt, e in rows:
    box(C0, t, W0, h, txt, GREY, GREYF); mids.append(mid(t, h)); t += h + gap
d_end = t - gap
box(C1, td, W1, d_end - td, 'no descriptor selection,\nso nothing to nest', GREY, GREYF, title='Control arms', center=True)
for ym in mids: ar(C0 + W0, ym, C1, ym)
cd = (td + d_end) / 2; hd3 = 13.0; sep = 9.5          # two models, each pool is modelled with both
tP1 = cd - sep - hd3 / 2; tR1 = cd + sep - hd3 / 2
box(C3, tP1, W3, hd3, 'PLS-DA\n1 latent variable\nleave-one-out', BLUE, BLUEF, bold=True)
box(C3, tR1, W3, hd3, 'Random forest\n300 trees\nout-of-bag', BLUE, BLUEF, bold=True)
assert tP1 >= td and tR1 + hd3 <= d_end
ar(C1 + W1, Y(cd - sep), C3, Y(cd - sep)); ar(C1 + W1, Y(cd + sep), C3, Y(cd + sep))
# audit round 2 (EW 04.10.2026): the Morgan2 arm WITH univariate filter is label-dependent and nested, so it is not a
# control arm but a sensitivity analysis -> separate orange branch below the control arms, own model (PLS-DA).
tF2, hF2 = d_end + 3.0, 18.0
box(C1, tF2, W1, hF2, 'nested like Route A,\nMorgan2 bits only\n(sensitivity analysis)', ORANGE, ORANGEF, title='Univariate filter')
yF2 = mid(tF2, hF2); xE7 = C0 + W0 / 2
line([xE7, xE7], [Y(d_end), yF2]); ar(xE7, yF2, C1, yF2)          # from the lower edge of the Morgan2 box
hd4 = 9.0
box(C3, tF2 + (hF2 - hd4) / 2, W3, hd4, 'PLS-DA', BLUE, BLUEF, bold=True)
ar(C1 + W1, yF2, C3, yF2)
d_end = tF2 + hF2

# ============================================ (d) replication and generalisation
te = d_end + 3.5; panel(te, 'd', 'Replication and generalisation'); te += 5.0
rowsd = [(12.0, '37 all-D peptides,\n4 strains,\n1,275 descriptors'), (12.0, 'Gaucher SELDI-TOF,\n39 sera,\n590 m/z variables'),
         (15.5, 'ChEMBL benchmark:\n16 series and\n4 size subsets,\nn = 25 to 60')]
gapd = 2.0; t = te; midsd = []
for h, txt in rowsd:
    box(C0, t, W0, h, txt, GREY, 'white'); midsd.append(mid(t, h)); t += h + gapd
de = t - gapd
box(C1, te, W1 + (C2 - C1 - W1) + W2, de - te, 'Routes A1, A1u and A2 with nested leave-one-out\nvalidation and whole-pipeline permutation.\nBenchmark: separate implementation\n(Methods, ChEMBL benchmark). Control arms\nfor the peptides and the benchmark.', GREY, GREYF, title='Same analysis plan', center=True)
for ym in midsd: ar(C0 + W0, ym, C1, ym)
# result boxes (EW 04:12 SAST): equal height, centres (and arrows) at 1/4 and 3/4 of the height of the box "Same pipeline",
# that is symmetric to its centre, equal margin at the top and at the bottom
hr = rowsd[2][0]; Hd = de - te
cr1, cr2 = te + Hd / 4, te + 3 * Hd / 4
tr1, tr2 = cr1 - hr / 2, cr2 - hr / 2
assert tr1 >= te and tr2 + hr <= de and tr1 + hr < tr2
box(C3, tr1, W3, hr, 'replication in\nindependent data', GREY, 'white')
box(C3, tr2, W3, hr, 'generalisation:\nfour pre-specified\nhypotheses H1 to H4', GREY, 'white')
xS = C2 + W2
ar(xS, Y(cr1), C3, Y(cr1))
ar(xS, Y(cr2), C3, Y(cr2))
NOTES.append(note((C1 + C2 + W2) / 2, de + 0.8, 'benchmark: selection rules, pipeline and hypotheses frozen before any model was fitted', GREY))
tE = de + 0.8 + 3.6; hE = 0.0

# ============================================ legend (2 rows of 3 entries)
tL = tE + hE + 3.0
LEG = [(GREY, 'white', 'data or result'), (GREEN, GREENF, 'label-free'), (ORANGE, GREENF, 'label-selected input'),
       (ORANGE, ORANGEF, 'label-dependent'), (BLUE, BLUEF, 'model'), (GREY, GREYF, 'control arm or further analysis')]
LEGT = []
for i, (ec, fc, lb) in enumerate(LEG):
    x0 = C0 + (i % 3) * 58.0; top = tL + (i // 3) * 5.6
    ax.add_patch(FancyBboxPatch((x0, Y(top) - 3.4), 4.6, 3.4, boxstyle='round,pad=0,rounding_size=0.8', linewidth=0.9,
                                edgecolor=ec, facecolor=fc))
    LEGT.append(ax.text(x0 + 6.0, Y(top) - 1.7, lb, fontsize=FN, va='center'))
L_END = tL + 5.6 + 3.4

# ============================================ checks
fig.canvas.draw(); r = fig.canvas.get_renderer()
bad, inner = [], []
for own, p in TXT:
    bp = p.get_window_extent(r)
    bbs = [t.get_window_extent(r) for t in own]
    for t, bt in zip(own, bbs):
        if not (bt.x0 >= bp.x0 + 2 and bt.x1 <= bp.x1 - 2 and bt.y0 >= bp.y0 + 1 and bt.y1 <= bp.y1 - 1):
            bad.append((t.get_text()[:30], round(bt.x0 - bp.x0), round(bp.x1 - bt.x1), round(bt.y0 - bp.y0), round(bp.y1 - bt.y1)))
    for i in range(len(bbs)):
        for j in range(i + 1, len(bbs)):
            if bbs[i].overlaps(bbs[j]): inner.append((own[i].get_text()[:20], own[j].get_text()[:20]))
ext = [p.get_window_extent(r) for p, _ in BOX]
ov = [(BOX[i][1], BOX[j][1]) for i in range(len(ext)) for j in range(i + 1, len(ext)) if ext[i].overlaps(ext[j])]
lov = [(L.get_text(), BOX[j][1]) for L in LETTERS for j, e in enumerate(ext) if L.get_window_extent(r).overlaps(e)]
nov = [(n.get_text()[:25], BOX[j][1]) for n in NOTES for j, e in enumerate(ext) if n.get_window_extent(r).overlaps(e)]
small = sorted({(round(t.get_fontsize(), 2), t.get_text()[:20]) for t in fig.findobj(matplotlib.text.Text)
                if t.get_text().strip() and t.get_fontsize() < MIN_PT})
figbb = fig.bbox
outside = [t.get_text()[:25] for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()
           and not (t.get_window_extent(r).x0 >= figbb.x0 and t.get_window_extent(r).x1 <= figbb.x1
                    and t.get_window_extent(r).y0 >= figbb.y0 and t.get_window_extent(r).y1 <= figbb.y1)]
legx = max(t.get_window_extent(r).x1 for t in LEGT) / r.points_to_pixels(72) * 25.4
assert not bad, ('text extends beyond the box', bad)
assert not inner, ('texts inside the box overlap', inner)
assert not ov, ('boxes overlap', ov)
assert not lov, ('panel letter overlaps a box', lov)
assert not nov, ('note overlaps a box', nov)
assert not small, ('font below 8 pt', small)
assert not outside, ('text outside the figure', outside)
assert L_END <= HMM - 0.5, ('height not sufficient', L_END, HMM)
assert abs(fig.get_size_inches()[0] * 25.4 - 174.0) < 0.01 and HMM <= 234.0
assert legx <= WMM - 1.0, legx
fonts = {t.get_fontname() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()}
print('Check passed | size %.0f x %.0f mm | content down to %.1f mm | smallest font %.1f pt | font %s | texts %d | boxes %d'
      % (WMM, HMM, L_END, min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()),
         fonts, sum(len(o) for o, _ in TXT), len(BOX)))

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else None
if OUTDIR:
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, 'Fig1_scheme_revised.png'), dpi=600, facecolor='white')
    fig.savefig(os.path.join(OUTDIR, 'Fig1_scheme_revised.pdf'), facecolor='white')
    print('Fig. 1 written to', OUTDIR)
    sys.exit(0)
if os.environ.get('FIG1_PREVIEW'):
    fig.savefig(os.environ['FIG1_PREVIEW'], dpi=300, facecolor='white'); sys.exit(0)
TS = datetime.datetime.now(zoneinfo.ZoneInfo('Europe/Berlin')).strftime('%d.%m.%Y %H%M')
for ext_ in ('png', 'pdf'):
    _osp.makedirs(_WORK + '/fig1', exist_ok=True)
    fig.savefig((_WORK + '/fig1/Fig1_Revision2_Entwurf %s.%s') % (TS, ext_), dpi=600 if ext_ == 'png' else None, facecolor='white')
print('Fig. 1 written', TS)
