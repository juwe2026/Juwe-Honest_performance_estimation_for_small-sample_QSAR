# -*- coding: utf-8 -*-
"""Fig. 1, Fassung 04.10.2026 (Wunsch EW): Druckgroesse nach den Springer-Vorgaben fuer Abbildungen.

Springer artwork guidelines: Schrift 8-12 pt in der gedruckten Groesse, Abbildungsbreite hoechstens 174 mm
(zweispaltig), Hoehe hoechstens 234 mm. Die Vorfassung (fig1_v2.py, durch diese Datei ersetzt) war 180 mm breit, mit 6-7 pt Schrift.
Diese Fassung zeichnet direkt in Millimetern auf einer 174 mm breiten Flaeche; jede Schrift ist >= 8 pt
(Text 8 pt, Kastentitel 8.5 pt fett, Feldbuchstaben 10 pt). Schriftart Liberation Sans (metrisch gleich Arial).
Inhalt, Farben, Pfeile und Aufbau wie fig1_v2.py; geaendert nur:
  - Zeilenumbrueche in den Kaesten (mehr Zeilen, weil die Schrift groesser ist),
  - E5/E6/E7 als Eckmarken in den Eingabekaesten von Panel d (wie E1-E4), nicht mehr auf den Pfeilen,
  - Legendenzeile zweizeilig (3 + 3 Felder),
  - die Gruppierung der Ueberlebenden unter Route A2 gruen mit orangem Rahmen ("input selected with labels").
  - Panel d: der Morgan2-Arm mit genestetem univariatem Filter als eigener oranger Zweig mit eigenem PLS-DA-Kasten
    (Sensitivitaetsanalyse, kein Kontrollarm; Audit Runde 2, 04.10.2026).
Fachlich: Route A2 filtert alle 2,208 Deskriptoren zuerst und gruppiert danach nur die Ueberlebenden; deshalb
fuehrt der Pfeil zu Route A2 vom Kasten der 2,208 Deskriptoren aus, und ihre Gruppierung liegt in der Schleife.
Pruefungen am Ende: jeder Text liegt in seinem Kasten, Texte eines Kastens ueberlappen nicht, Kaesten und
Feldbuchstaben ueberlappen nicht, jede Schrift >= 8 pt, Breite = 174 mm, Hoehe <= 234 mm.

Usage in the code archive:  cd figures && python ../revision/fig1_v3.py output
  writes output/Fig1_scheme_revised.png (600 dpi, as printed in the revised article) and output/Fig1_scheme_revised.pdf.
Without an argument the script writes the time-stamped working files of the revision instead."""
import sys, os, datetime, zoneinfo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = 'Liberation Sans'
plt.rcParams['pdf.fonttype'] = 42          # TrueType eingebettet
plt.rcParams['svg.fonttype'] = 'path'

GREEN = '#2E7D32'; GREENF = '#E7F3E8'; ORANGE = '#D96A00'; ORANGEF = '#FCECD9'
BLUE = '#1A5FB4'; BLUEF = '#E3EBF8'; GREY = '#5A5A5A'; GREYF = '#EFEFEF'; LINE = '#505050'
FB, FT, FTAG, FN, FP = 8.0, 8.5, 8.0, 8.0, 10.0      # Text, Kastentitel, E-Marken, Hinweise, Feldbuchstaben (pt)
MIN_PT = 8.0
LS = 1.25                                             # Zeilenabstand
WMM, HMM = 174.0, 226.0                               # Abbildungsgroesse in mm
PAD = 1.8                                             # Innenabstand oben (mm)
LT = 4.3                                              # Hoehe der Titelzeile (mm)

fig = plt.figure(figsize=(WMM / 25.4, HMM / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, WMM); ax.set_ylim(0, HMM); ax.axis('off')
def Y(top): return HMM - top                          # Abstand von oben (mm) -> y-Koordinate

# Spalten (mm): Luecke nach Spalte 0 breiter wegen der Winkelpfeile
C0, W0 = 2.0, 35.0
C1, W1 = 48.0, 38.0
C2, W2 = 93.0, 38.0
C3, W3 = 138.0, 34.0
BOX = []; TXT = []; LETTERS = []

def box(x, top, w, h, txt, ec, fc, bold=False, title=None, tag=None, tc=None):
    """Kasten mit oberer Kante 'top' (mm von oben). tc: Titelfarbe, falls sie vom Rahmen abweicht."""
    y = Y(top) - h
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.5', linewidth=0.9,
                       edgecolor=ec, facecolor=fc, zorder=3)
    ax.add_patch(p); BOX.append((p, txt[:20]))
    own = []
    if title:
        own.append(ax.text(x + w / 2, Y(top) - PAD, title, ha='center', va='top', fontsize=FT, fontweight='bold',
                           color=tc or ec, zorder=4))
        own.append(ax.text(x + w / 2, Y(top) - PAD - LT, txt, ha='center', va='top', fontsize=FB, color='#1A1A1A',
                           zorder=4, linespacing=LS))
    else:
        own.append(ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=FB, color='#1A1A1A', zorder=4,
                           fontweight='bold' if bold else 'normal', linespacing=LS))
    if tag:
        own.append(ax.text(x + w - 1.0, Y(top) - 0.9, tag, ha='right', va='top', fontsize=FTAG, color='#666666',
                           fontweight='bold', zorder=5))
    TXT.append((own, p))
    return p

def ar(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=8, linewidth=0.85, color=LINE,
                                 zorder=2, shrinkA=0, shrinkB=0))

def line(xs, ys):
    ax.add_line(Line2D(xs, ys, linewidth=0.85, color=LINE, zorder=2, solid_capstyle='butt', solid_joinstyle='miter'))

def elbow(x1, y1, xm, y2, x2):
    """waagrecht nach rechts bis xm, senkrecht auf y2, dann waagrecht mit Pfeil bis x2"""
    line([x1, xm], [y1, y1]); line([xm, xm], [y1, y2]); ar(xm, y2, x2, y2)

def panel(top, t):
    LETTERS.append(ax.text(0.2, Y(top), t, fontsize=FP, fontweight='bold', color='#1A1A1A', ha='left', va='top', zorder=6))

def note(xc, top, t, color):
    return ax.text(xc, Y(top), t, ha='center', va='top', fontsize=FN, color=color, style='italic', linespacing=LS)

def mid(top, h): return Y(top + h / 2)

# ============================================ (a) Route A
panel(0.5, 'a')
FX0, FX1 = C1 - 1.8, C2 + W2 + 1.8
FT_, FB_ = 11.5, 69.5                                  # Rahmen oben/unten (mm von oben)
ax.add_patch(Rectangle((FX0, Y(FB_)), FX1 - FX0, FB_ - FT_, linewidth=1.0, edgecolor=ORANGE, facecolor='none',
                       linestyle=(0, (4, 2.6)), zorder=1))
NOTES = [ax.text((FX0 + FX1) / 2, Y(FT_) + 1.0, 'repeated inside every leave-one-out fold and\nunder each of the 2,000 label permutations',
                 ha='center', va='bottom', fontsize=FN, color=ORANGE, style='italic', linespacing=LS)]
tD, hD = 13.0, 9.5;  tP, hP = 27.0, 13.5;  tG, hG = 45.0, 21.0
box(C0, tD, W0, hD, '31 monoterpenoids,\n5 binary endpoints', GREY, 'white', bold=True)
box(C0, tP, W0, hP, '2,208 cleaned\ndescriptors from\n5 libraries', GREEN, GREENF)
box(C0, tG, W0, hG, '|r| or |rho| > 0.7\n2,208 → 182 groups\nties → largest variance', GREEN, GREENF, title='Correlation grouping')
ar(C0 + W0 / 2, Y(tD + hD), C0 + W0 / 2, Y(tP)); ar(C0 + W0 / 2, Y(tP + hP), C0 + W0 / 2, Y(tG))
NOTES.append(note(C0 + W0 / 2, tG + hG + 0.8, 'label-free, computed once', GREEN))
tA2, hA2 = 13.0, 18.0;  tA1, hA1 = 34.5, 14.5;  tA1u, hA1u = 52.5, 14.5
box(C1, tA2, W1, hA2, 'univariate filter\napplied to all\n2,208 descriptors', ORANGE, ORANGEF, title='Route A2', tag='E1')
box(C1, tA1, W1, hA1, 'representative per group\nchosen by effect size', ORANGE, ORANGEF, title='Route A1', tag='E1')
box(C1, tA1u, W1, hA1u, 'representative per group\nchosen by centrality', GREEN, GREENF, title='Route A1u', tag='E1')
yA2, yA1, yA1u = mid(tA2, hA2), mid(tA1, hA1), mid(tA1u, hA1u)
assert Y(tG + hG) < yA1u < Y(tG), 'A1u-Pfeil muss aus dem Gruppierungskasten kommen'
# Pfeile ohne Kreuzung: 2,208 Deskriptoren -> Route A2 (Filter zuerst); Gruppierung -> Route A1 und A1u
xA2m, xGm = C0 + W0 + 3.0, C0 + W0 + 6.5
elbow(C0 + W0, mid(tP, hP), xA2m, yA2, C1)
ar(C0 + W0, yA1u, C1, yA1u)
line([xGm, xGm], [yA1u, yA1]); ar(xGm, yA1, C1, yA1)
ax.plot([xGm], [yA1u], marker='o', markersize=2.2, color=LINE, zorder=3)
# Gruppierung der Ueberlebenden: Schritt label-frei (gruen), Eingabe aus dem label-gestuetzten Filter (oranger Rahmen)
tS, hS = 13.0, 18.0;  tF, hF = 38.0, 18.0
box(C2, tS, W2, hS, 'of the survivors;\ninput selected\nwith labels', ORANGE, GREENF, title='Correlation grouping', tc=GREEN)
box(C2, tF, W2, hF, 'Mann–Whitney p ≤ 0.05;\nthe |r_rb| > 0.3\nthreshold never binds', ORANGE, ORANGEF, title='Univariate filter')
assert Y(tS + hS) < yA2 < Y(tS)
ar(C1 + W1, yA2, C2, yA2)
ar(C1 + W1, yA1, C2, Y(tF + 5.0))
ar(C1 + W1, yA1u, C2, Y(tF + hF - 4.0))
tM, hM = 18.0, 14.5
box(C3, tM, W3, hM, 'PLS-DA\n1 latent variable\nleave-one-out', BLUE, BLUEF, bold=True)
assert Y(tM + hM) < yA2 < Y(tM)
ar(C2 + W2, yA2, C3, yA2)                              # Gruppierung der Ueberlebenden -> Modell (waagrecht)
# symmetrisch zur Mitte des Filterkastens: Start +/- 2 mm, Ende +/- d
yFc = mid(tF, hF); d = Y(tM + hM - 2.5) - yFc       # > 0: Ende oben im Modellkasten, unten symmetrisch
hN = 11.0; tN = (HMM - (yFc - d)) - hN / 2; assert d > 0
box(C3, tN, W3, hN, 'Null distribution\nof the hit count', ORANGE, ORANGEF, tag='E4')
ar(C2 + W2, yFc + 2.0, C3, yFc + d)
ar(C2 + W2, yFc - 2.0, C3, yFc - d)
A_END = max(FB_, tN + hN, tG + hG + 5.0)

# ============================================ (b) Route B
tb = A_END + 3.5; panel(tb, 'b'); tb += 5.0
hb0, hb1, hb3 = 13.5, 14.5, 10.0; cb = tb + hb1 / 2
box(C0, cb - hb0 / 2, W0, hb0, '25 literature\ndescriptors\nreduced to 12', GREEN, GREENF)
box(C1, cb - hb1 / 2, W1, hb1, 'no supervised selection,\nso nothing to nest', GREEN, GREENF, title='Route B', tag='E2')
box(C3, cb - hb3 / 2, W3, hb3, 'PLS-DA\n1 and 3 LV', BLUE, BLUEF, bold=True)
ar(C0 + W0, Y(cb), C1, Y(cb)); ar(C1 + W1, Y(cb), C3, Y(cb))
B_END = cb + hb1 / 2

# ============================================ (c) two-stage VIP
tc_ = B_END + 3.5; panel(tc_, 'c'); tc_ += 5.0
hc0, hc1, hc3 = 13.5, 14.5, 10.0; cc = tc_ + hc1 / 2
box(C0, cc - hc0 / 2, W0, hc0, 'Route A1\ndescriptor set\n(35 to 56)', ORANGE, ORANGEF)
box(C1, cc - hc1 / 2, W1, hc1, 'ranks the descriptors\nby VIP score', BLUE, BLUEF, title='First PLS-DA', tag='E3')
box(C2, cc - hc1 / 2, W2, hc1, 'on the 15 highest-ranked\ndescriptors', BLUE, BLUEF, title='Second PLS-DA')
box(C3, cc - hc3 / 2, W3, hc3, 'reported\nperformance', GREY, 'white')
for xa, xb in ((C0 + W0, C1), (C1 + W1, C2), (C2 + W2, C3)): ar(xa, Y(cc), xb, Y(cc))
NOTES.append(note((C1 + C3 + W3) / 2, cc + hc1 / 2 + 0.8,
                  'label-dependent throughout; nested and apparent variants\nare compared in Tables 3 and 4', ORANGE))
C_END = cc + hc1 / 2 + 0.8 + 2 * 3.6

# ============================================ (d) control arms
td = C_END + 3.0; panel(td, 'd'); td += 5.0
rows = [(10.5, 'all 2,208\ndescriptors', 'E5'), (14.0, '182 label-free\nrepresentatives\n(panel a)', 'E6'),
        (10.5, 'Morgan2\nfingerprint bits', 'E7')]
gap = 2.5; t = td; mids = []
for h, txt, e in rows:
    box(C0, t, W0, h, txt, GREY, GREYF, tag=e); mids.append(mid(t, h)); t += h + gap
d_end = t - gap
box(C1, td, W1, d_end - td, 'no descriptor selection,\nso nothing to nest', GREY, GREYF, title='Control arms')
for ym in mids: ar(C0 + W0, ym, C1, ym)
cd = (td + d_end) / 2; hd3 = 14.0
box(C3, cd - hd3 / 2, W3, hd3, 'PLS-DA and\nrandom forest\n(out-of-bag)', BLUE, BLUEF)
ar(C1 + W1, Y(cd), C3, Y(cd))
# Audit Runde 2 (EW 04.10.2026): der Morgan2-Arm MIT univariatem Filter ist label-abhaengig und genestet, also kein
# Kontrollarm, sondern eine Sensitivitaetsanalyse -> eigener oranger Zweig unter den Kontrollarmen, eigenes Modell (PLS-DA).
tF2, hF2 = d_end + 3.0, 18.0
box(C1, tF2, W1, hF2, 'nested like Route A;\nMorgan2 bits only\n(sensitivity analysis)', ORANGE, ORANGEF, title='Univariate filter')
yF2 = mid(tF2, hF2); xE7 = C0 + W0 / 2
line([xE7, xE7], [Y(d_end), yF2]); ar(xE7, yF2, C1, yF2)          # von der Unterkante des Morgan2-Kastens
hd4 = 9.0
box(C3, tF2 + (hF2 - hd4) / 2, W3, hd4, 'PLS-DA', BLUE, BLUEF, bold=True)
ar(C1 + W1, yF2, C3, yF2)
d_end = tF2 + hF2

# ============================================ E8, E9
tE = d_end + 4.0; hE = 14.5
WH = (C3 + W3 - C0 - 2.5) / 2
box(C0, tE, WH, hE, 'E8  Sensitivity analyses: number of latent variables,\ntie-break rule, decision threshold, 3D descriptors,\n'
    'fold-internal grouping and scaling, resampling', GREY, GREYF)
box(C0 + WH + 2.5, tE, WH, hE, 'E9  Identical pipeline on two independent data sets:\nGaucher SELDI-TOF (590 m/z) and 37 all-D\n'
    'octapeptides (1,275 descriptors)', GREY, GREYF)

# ============================================ Legende (2 Zeilen zu 3 Feldern)
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

# ============================================ Pruefungen
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
assert not bad, ('Text ragt aus dem Kasten', bad)
assert not inner, ('Texte im Kasten ueberlappen', inner)
assert not ov, ('Kaesten ueberlappen', ov)
assert not lov, ('Feldbuchstabe ueberlappt Kasten', lov)
assert not nov, ('Hinweistext ueberlappt Kasten', nov)
assert not small, ('Schrift unter 8 pt', small)
assert not outside, ('Text ausserhalb der Abbildung', outside)
assert L_END <= HMM - 0.5, ('Hoehe reicht nicht', L_END, HMM)
assert abs(fig.get_size_inches()[0] * 25.4 - 174.0) < 0.01 and HMM <= 234.0
assert legx <= WMM - 1.0, legx
fonts = {t.get_fontname() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()}
print('Pruefung bestanden | Groesse %.0f x %.0f mm | Inhalt bis %.1f mm | kleinste Schrift %.1f pt | Schrift %s | Texte %d | Kaesten %d'
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
for ext_ in ('png', 'pdf', 'svg'):
    fig.savefig('/mnt/user-data/outputs/Fig1_neu_mit_Panels %s.%s' % (TS, ext_), dpi=600 if ext_ == 'png' else None, facecolor='white')
fig.savefig('/home/claude/rev/figs_rev/Fig1_preview.png', dpi=200, facecolor='white')
print('Fig. 1 geschrieben', TS)
