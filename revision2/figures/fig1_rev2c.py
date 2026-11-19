# -*- coding: utf-8 -*-
"""
Change of 17.11.2026 (approved by EW, item 6 of 17.11.2026): the box of the univariate filter reads
"the |effect size| > 0.3" instead of "the |r_rb| > 0.3", whose subscript was printed as a literal underscore. The
effect size is the rank-biserial correlation defined in the Methods. The approved wording in fig1_rev2.py was changed
in the same way. Nothing else changed.

Change of 14.11.2026 (approved by EW, item V8a): the dashed frame no longer encloses the box of Route A1u. The
representative chosen by centrality is computed once on all compounds without the labels (unsup_rep.centrality_reps),
so it is not one of the selection steps repeated inside every fold. Only the outline of the frame changed. All boxes,
arrows, texts and the size of the figure are unchanged.
Change of 07.11.2026 (approved by EW, items Ä9, Ä27 and F01 of the review of 06.11.2026): panel d names the
separate implementation of the benchmark (Methods, ChEMBL benchmark) and "37 all-D peptides" (six of the
37 peptides are shorter than eight residues), the note above the dashed frame names the selection steps.
Change of 07.11.2026 (approved by EW, variant A2): 186.5 mm instead of 199.5 mm high, so that figure and legend fit into
225 mm. Smaller gaps between boxes and panels, the note above the dashed frame in one line. Wording and font sizes
unchanged.
Fig. 1, Revision 2, version for the BMC requirements, with panels c and d one below the other as in the approved
version (EW 06.11.2026 18:09 Berlin). Panels a and b and the colour key as in the variant with c and d side by side
(not archived). In c and d the column of the inputs is wider (47 mm), so the boxes need fewer lines, and the models
stand in the model column of panels a and b.

Requirements (Journal of Cheminformatics, Figure size and resolution): full width 170 mm, at most 225 mm height for
figure and legend together, text legible at final size, lines > 0.25 pt.
Aim: 170 mm wide and every font >= 8 pt at this final size.

Content, texts, colours and arrows as the approved version fig1_rev2.py (06.11.2026 03:12). Changed only:
  - width 170 instead of 174 mm, columns of panels a and b correspondingly narrower,
  - in panels c and d a wider input column, models in the model column of panels a and b,
  - the branch of the sensitivity analysis (univariate filter on Morgan2 bits) as a box below the control arms,
    the line starts at the lower edge of the Morgan2 box, as approved,
  - line breaks in the boxes where the boxes are narrower or wider (wording unchanged, checked at the end).
Checks at the end as in fig1_rev2.py, in addition: the wording of all texts equals the approved version (white space and
line breaks excepted), width 170 mm.

Usage in the code archive:  cd figures && python ../revision2/figures/fig1_rev2c.py output
  writes output/Fig1_scheme_revised.png (600 dpi) and output/Fig1_scheme_revised.pdf.
Without an argument the script writes the time-stamped working files of the revision instead."""
import sys, os, re, datetime, zoneinfo, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = 'Liberation Sans'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'path'

GREEN = '#2E7D32'; GREENF = '#E7F3E8'; ORANGE = '#D96A00'; ORANGEF = '#FCECD9'
BLUE = '#1A5FB4'; BLUEF = '#E3EBF8'; GREY = '#5A5A5A'; GREYF = '#EFEFEF'; LINE = '#505050'
FB, FT, FTAG, FN, FP = 8.0, 8.5, 8.0, 8.0, 10.0
MIN_PT = 8.0
LS = 1.25
WMM = 170.0
PAD = 1.8
LT = 4.3
LH = FB * LS * 25.4 / 72                 # line height (mm)


def build(HMM):
    fig = plt.figure(figsize=(WMM / 25.4, HMM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, WMM); ax.set_ylim(0, HMM); ax.axis('off')
    Y = lambda top: HMM - top
    BOX, TXT, LETTERS, NOTES = [], [], [], []

    def box(x, top, w, h, txt, ec, fc, bold=False, title=None, tc=None, center=False):
        y = Y(top) - h
        p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.5', linewidth=0.9,
                           edgecolor=ec, facecolor=fc, zorder=3)
        ax.add_patch(p); BOX.append((p, txt[:20]))
        own = []
        if title:
            off = 0.0
            if center:
                blk = LT + txt.count('\n') * LH + FB * 25.4 / 72
                off = max(0.0, (h - 2 * PAD - blk) / 2)
            t0 = top + off
            own.append(ax.text(x + w / 2, Y(t0) - PAD, title, ha='center', va='top', fontsize=FT, fontweight='bold',
                               color=tc or ec, zorder=4))
            own.append(ax.text(x + w / 2, Y(t0) - PAD - LT, txt, ha='center', va='top', fontsize=FB, color='#1A1A1A',
                               zorder=4, linespacing=LS))
        else:
            own.append(ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=FB, color='#1A1A1A', zorder=4,
                               fontweight='bold' if bold else 'normal', linespacing=LS))
        TXT.append((own, p))
        return p

    def ar(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=8, linewidth=0.85, color=LINE,
                                     zorder=2, shrinkA=0, shrinkB=0))

    def line(xs, ys):
        ax.add_line(Line2D(xs, ys, linewidth=0.85, color=LINE, zorder=2, solid_capstyle='butt', solid_joinstyle='miter'))

    def elbow(x1, y1, xm, y2, x2):
        line([x1, xm], [y1, y1]); line([xm, xm], [y1, y2]); ar(xm, y2, x2, y2)

    def panel(x0, top, t, title=''):
        LETTERS.append(ax.text(x0 + 0.2, Y(top), t, fontsize=FP, fontweight='bold', color='#1A1A1A', ha='left', va='top', zorder=6))
        if title:
            LETTERS.append(ax.text(x0 + 4.2, Y(top) - 0.35, title, fontsize=FT, fontweight='bold', color='#1A1A1A',
                                   ha='left', va='top', zorder=6))

    def note(xc, top, t, color):
        n = ax.text(xc, Y(top), t, ha='center', va='top', fontsize=FN, color=color, style='italic', linespacing=LS)
        NOTES.append(n); return n

    def mid(top, h): return Y(top + h / 2)

    # ------------------------------------------------ columns of panels a and b (mm)
    C0, W0 = 1.5, 33.0
    C1, W1 = 45.0, 36.0
    C2, W2 = 88.5, 36.0
    C3, W3 = 132.0, 36.5

    # ============================================ (a) Route A
    panel(0, 0.5, 'a', 'Route A, worked example')
    FX0, FX1 = C1 - 1.8, C2 + W2 + 1.8
    FT_, FB_ = 8.5, 64.0                                  # FB_ only keeps the layout below panel a unchanged
    # the frame encloses Route A2, Route A1, the grouping of the survivors and the univariate filter, not Route A1u
    FM_, FXS, FB2 = 46.25, C2 - 3.0, 53.3                 # step between A1 and A1u, left edge and lower edge at the filter
    ax.add_patch(Polygon([(FX0, Y(FT_)), (FX1, Y(FT_)), (FX1, Y(FB2)), (FXS, Y(FB2)), (FXS, Y(FM_)), (FX0, Y(FM_))], closed=True,
                         linewidth=1.0, edgecolor=ORANGE, facecolor='none', linestyle=(0, (4, 2.6)), joinstyle='miter', zorder=1))
    NOTES.append(ax.text((FX0 + FX1) / 2, Y(FT_) + 1.0, 'selection steps repeated inside every leave-one-out fold and under each of the 2,000 label permutations',
                         ha='center', va='bottom', fontsize=FN, color=ORANGE, style='italic', linespacing=LS))
    tD, hD = 10.0, 9.5;  tP, hP = 23.0, 13.5;  tG, hG = 40.0, 19.5
    box(C0, tD, W0, hD, '31 monoterpenoids,\n5 binary endpoints', GREY, 'white', bold=True)
    box(C0, tP, W0, hP, '2,208 cleaned\ndescriptors from\n5 libraries', GREEN, GREENF)
    box(C0, tG, W0, hG, '|r| or |rho| > 0.7\n2,208 → 182 groups\nties → largest variance', GREEN, GREENF, title='Correlation grouping')
    ar(C0 + W0 / 2, Y(tD + hD), C0 + W0 / 2, Y(tP)); ar(C0 + W0 / 2, Y(tP + hP), C0 + W0 / 2, Y(tG))
    note(C0 + W0 / 2, tG + hG + 0.8, 'label-free, computed once', GREEN)
    tA2, hA2 = 10.0, 18.0;  tA1, hA1 = 30.5, 14.5;  tA1u, hA1u = 47.5, 14.5
    box(C1, tA2, W1, hA2, 'univariate filter\napplied to all\n2,208 descriptors', ORANGE, ORANGEF, title='Route A2')
    box(C1, tA1, W1, hA1, 'representative per group\nchosen by effect size', ORANGE, ORANGEF, title='Route A1')
    box(C1, tA1u, W1, hA1u, 'representative per group\nchosen by centrality', GREEN, GREENF, title='Route A1u')
    yA2, yA1, yA1u = mid(tA2, hA2), mid(tA1, hA1), mid(tA1u, hA1u)
    assert tA1 + hA1 + 1.0 <= FM_ <= tA1u - 1.0 and C1 + W1 + 1.0 <= FXS <= C2 - 1.0      # the frame clears A1, A1u and the filter
    assert Y(tG + hG) < yA1u < Y(tG)
    xA2m, xGm = C0 + W0 + 3.0, C0 + W0 + 6.5
    elbow(C0 + W0, mid(tP, hP), xA2m, yA2, C1)
    ar(C0 + W0, yA1u, C1, yA1u)
    line([xGm, xGm], [yA1u, yA1]); ar(xGm, yA1, C1, yA1)
    ax.plot([xGm], [yA1u], marker='o', markersize=2.2, color=LINE, zorder=3)
    tS, hS = 10.0, 18.0;  tF, hF = 33.5, 18.0
    box(C2, tS, W2, hS, 'of the survivors,\ninput selected\nwith labels', ORANGE, GREENF, title='Correlation grouping', tc=GREEN)
    assert tF + hF + 1.0 <= FB2 and FB2 < tA1u + hA1u
    box(C2, tF, W2, hF, 'Mann–Whitney p ≤ 0.05,\nthe |effect size| > 0.3\nthreshold never binds', ORANGE, ORANGEF, title='Univariate filter')
    assert Y(tS + hS) < yA2 < Y(tS)
    ar(C1 + W1, yA2, C2, yA2)
    ar(C1 + W1, yA1, C2, Y(tF + 5.0))
    ar(C1 + W1, yA1u, C2, Y(tF + hF - 4.0))
    tM, hM = 15.0, 14.5
    box(C3, tM, W3, hM, 'PLS-DA\n1 latent variable\nleave-one-out', BLUE, BLUEF, bold=True)
    assert Y(tM + hM) < yA2 < Y(tM)
    ar(C2 + W2, yA2, C3, yA2)
    yFc = mid(tF, hF); d = Y(tM + hM - 2.5) - yFc
    hN = 11.0; tN = (HMM - (yFc - d)) - hN / 2; assert d > 0
    box(C3, tN, W3, hN, 'Null distribution\nof the hit count', ORANGE, ORANGEF)
    ar(C2 + W2, yFc + 2.0, C3, yFc + d)
    ar(C2 + W2, yFc - 2.0, C3, yFc - d)
    A_END = max(FB_, tN + hN, tG + hG + 5.0)

    # ============================================ (b) Route B
    tb = A_END + 2.5; panel(0, tb, 'b', 'Route B'); tb += 4.5
    hb0, hb1, hb3 = 13.5, 14.5, 10.0; cb = tb + hb1 / 2
    box(C0, cb - hb0 / 2, W0, hb0, '25 literature\ndescriptors\nreduced to 12', GREEN, GREENF)
    box(C1, cb - hb1 / 2, W1, hb1, 'no supervised selection,\nso nothing to nest', GREEN, GREENF, title='Route B')
    box(C3, cb - hb3 / 2, W3, hb3, 'PLS-DA\n1 and 3 LV', BLUE, BLUEF, bold=True)
    ar(C0 + W0, Y(cb), C1, Y(cb)); ar(C1 + W1, Y(cb), C3, Y(cb))
    B_END = cb + hb1 / 2

    # ============================================ (c) control arms, full width, below panel b (as approved)
    XI, WI = 1.5, 47.0                                    # input column of c and d (wider than C0, so fewer lines)
    XK, WK = 58.0, 36.0                                   # control arms box
    tc = B_END + 2.0; panel(0, tc, 'c', 'Control arms and a sensitivity branch'); tc += 4.5
    rows = [(6.5, 'all 2,208 descriptors'), (10.0, '182 label-free representatives\n(panel a)'),
            (6.5, 'Morgan2 fingerprint bits')]
    gap = 1.0; t = tc; mids = []
    for h, txt in rows:
        box(XI, t, WI, h, txt, GREY, GREYF); mids.append(mid(t, h)); t += h + gap
    c_in_end = t - gap
    box(XK, tc, WK, c_in_end - tc, 'no descriptor selection,\nso nothing to nest', GREY, GREYF, title='Control arms', center=True)
    for ym in mids: ar(XI + WI, ym, XK, ym)
    cdm = (tc + c_in_end) / 2; hd3 = 11.5; sep = 6.75    # two models, each pool is modelled with both
    tP1 = cdm - sep - hd3 / 2; tR1 = cdm + sep - hd3 / 2
    box(C3, tP1, W3, hd3, 'PLS-DA\n1 latent variable\nleave-one-out', BLUE, BLUEF, bold=True)
    box(C3, tR1, W3, hd3, 'Random forest\n300 trees\nout-of-bag', BLUE, BLUEF, bold=True)
    assert tP1 >= tc and tR1 + hd3 <= c_in_end and tP1 + hd3 < tR1, (tP1, tc, tR1 + hd3, c_in_end)
    ar(XK + WK, Y(cdm - sep), C3, Y(cdm - sep)); ar(XK + WK, Y(cdm + sep), C3, Y(cdm + sep))
    # sensitivity branch (Morgan2 bits with the nested univariate filter), as approved: from the bottom of the Morgan2 box
    tF2, hF2 = c_in_end + 2.0, 14.5
    WF2 = 52.0
    box(XK, tF2, WF2, hF2, 'nested like Route A, Morgan2 bits only\n(sensitivity analysis)', ORANGE, ORANGEF, title='Univariate filter')
    yF2 = mid(tF2, hF2); xM = XI + WI / 2
    line([xM, xM], [Y(c_in_end), yF2]); ar(xM, yF2, XK, yF2)
    hd4 = 9.0
    box(C3, tF2 + (hF2 - hd4) / 2, W3, hd4, 'PLS-DA', BLUE, BLUEF, bold=True)
    ar(XK + WF2, yF2, C3, yF2)
    C_END = tF2 + hF2

    # ============================================ (d) replication and generalisation, full width, below panel c
    td = C_END + 2.0; panel(0, td, 'd', 'Replication and generalisation'); td += 4.5
    rowsd = [(9.5, '37 all-D peptides, 4 strains,\n1,275 descriptors'), (9.5, 'Gaucher SELDI-TOF, 39 sera,\n590 m/z variables'),
             (9.5, 'ChEMBL benchmark: 16 series and\n4 size subsets, n = 25 to 60')]
    gapd = 1.0; t = td; midsd = []
    for h, txt in rowsd:
        box(XI, t, WI, h, txt, GREY, 'white'); midsd.append(mid(t, h)); t += h + gapd
    de = t - gapd
    XS, WS = XK, C2 + W2 - XK
    box(XS, td, WS, de - td, 'Routes A1, A1u and A2 with nested leave-one-out\nvalidation and whole-pipeline permutation.\nBenchmark: separate implementation\n(Methods, ChEMBL benchmark). Control arms\nfor the peptides and the benchmark.',
        GREY, GREYF, title='Same analysis plan', center=True)
    for ym in midsd: ar(XI + WI, ym, XS, ym)
    # result boxes as approved: equal height, centres (and arrows) at 1/4 and 3/4 of the height of "Same pipeline"
    hr = 13.0; Hd = de - td
    cr1, cr2 = td + Hd / 4, td + 3 * Hd / 4
    tr1, tr2 = cr1 - hr / 2, cr2 - hr / 2
    assert tr1 >= td and tr2 + hr <= de and tr1 + hr < tr2, (tr1, td, tr2 + hr, de)
    box(C3, tr1, W3, hr, 'replication in\nindependent data', GREY, 'white')
    box(C3, tr2, W3, hr, 'generalisation:\nfour pre-specified\nhypotheses H1 to H4', GREY, 'white')
    ar(XS + WS, Y(cr1), C3, Y(cr1)); ar(XS + WS, Y(cr2), C3, Y(cr2))
    note(XS + WS / 2, de + 0.8, 'benchmark: selection rules, pipeline and hypotheses frozen before any model was fitted', GREY)
    D_END = de + 0.8 + FN * 25.4 / 72 * 1.2 + 0.4
    C_END = D_END

    # ============================================ legend (2 rows of 3)
    tL = D_END + 2.0
    LEG = [(GREY, 'white', 'data or result'), (GREEN, GREENF, 'label-free'), (ORANGE, GREENF, 'label-selected input'),
           (ORANGE, ORANGEF, 'label-dependent'), (BLUE, BLUEF, 'model'), (GREY, GREYF, 'control arm or further analysis')]
    LEGT = []
    for i, (ec, fc, lb) in enumerate(LEG):
        x0 = C0 + (i % 3) * 57.0; top = tL + (i // 3) * 4.6
        ax.add_patch(FancyBboxPatch((x0, Y(top) - 3.4), 4.6, 3.4, boxstyle='round,pad=0,rounding_size=0.8', linewidth=0.9,
                                    edgecolor=ec, facecolor=fc))
        LEGT.append(ax.text(x0 + 6.0, Y(top) - 1.7, lb, fontsize=FN, va='center'))
    L_END = tL + 4.6 + 3.4
    return fig, ax, BOX, TXT, LETTERS, NOTES, LEGT, L_END


# first pass: height of the content, second pass: figure exactly that high
fig, *_, L_END = build(250.0); plt.close(fig)
HMM = math.ceil((L_END + 0.6) * 2) / 2
fig, ax, BOX, TXT, LETTERS, NOTES, LEGT, L_END = build(HMM)

# ============================================ checks (as fig1_rev2.py)
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
assert abs(fig.get_size_inches()[0] * 25.4 - WMM) < 0.01
assert legx <= WMM - 1.0, legx
# wording equal to the approved version (white space and line breaks excepted)
APPROVED = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig1_rev2.py')
if os.path.exists(APPROVED):
    src = open(APPROVED, encoding='utf-8').read()
    src = src[src.index('# ============================================ (a) Route A'):src.index('# ============================================ checks')]
    src = '\n'.join(l for l in src.split('\n') if not l.lstrip().startswith('assert') and not l.lstrip().startswith('#'))
    norm = lambda s: re.sub(r'\s+', ' ', s.replace('\\n', ' ')).strip()
    want = sorted(norm(m) for m in re.findall(r"'([^'\n]*?(?:\\n[^'\n]*?)*)'", src) if re.search(r'[a-z]{3}', m) and ' ' in m)
    have = sorted(norm(t.get_text()) for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and ' ' in t.get_text())
    miss = [w for w in want if w not in have]; extra = [h for h in have if h not in want]
    assert not miss and not extra, ('wording differs', miss, extra)
    print('wording equal to the approved version:', len(have), 'texts')
fonts = {t.get_fontname() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()}
print('Check passed | size %.0f x %.1f mm | content down to %.1f mm | smallest font %.1f pt | font %s | texts %d | boxes %d'
      % (WMM, HMM, L_END, min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()),
         fonts, sum(len(o) for o, _ in TXT), len(BOX)))

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else None
if OUTDIR:
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, 'Fig1_scheme_revised.png'), dpi=600, facecolor='white')
    fig.savefig(os.path.join(OUTDIR, 'Fig1_scheme_revised.pdf'), facecolor='white')
    print('Fig. 1 written to', OUTDIR); sys.exit(0)
if os.environ.get('FIG1_PREVIEW'):
    fig.savefig(os.environ['FIG1_PREVIEW'], dpi=200, facecolor='white'); sys.exit(0)
TS = datetime.datetime.now(zoneinfo.ZoneInfo('Europe/Berlin')).strftime('%d.%m.%Y %H%M')
for ext_ in ('png', 'pdf'):
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Fig1_Revision2_BMC_cd_untereinander_Entwurf %s.%s' % (TS, ext_)),
                dpi=600 if ext_ == 'png' else None, facecolor='white')
print('Fig. 1 written', TS)
