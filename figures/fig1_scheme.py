"""Fig. 1: the descriptor-selection routes evaluated in this study.

Neu aufgebaut am 22. August, ueberarbeitet am 23. August nach der Autorenpruefung
des Word-Ausdrucks. Zwei Probleme wurden behoben:

1. Cleaning und Descriptor pool hatten unterschiedliche Kastenhoehen, sodass die
   Verteilerschiene zu den drei Routen durch den Kasten Descriptor pool lief.
   Beide Kaesten haben jetzt dieselbe Hoehe.
2. Final set, PLS-DA, Nested validation required (Route A), Simple permutation
   is already honest (Route B) und der VIP-Hinweis (Route C) standen links bzw.
   unterhalb der Ketten und liessen viel Weissraum rechts frei. Sie stehen jetzt
   am rechten Rand der Abbildung, Final set/PLS-DA auf Hoehe von Route A1u,
   Nested validation required unterhalb von Route A2, und die drei Hinweiskaesten
   in Route A/B/C haben dieselbe Breite.

Die Breite jedes Kastens wird weiterhin aus der laengsten Textzeile berechnet, und
die Gesamtabbildung wird am Ende auf die tatsaechlich benutzte Breite und Hoehe
zugeschnitten (nicht auf eine angenommene Breite von 100 Einheiten), damit nichts
abgeschnitten wird und die Schrift durchgehend im richtigen Massstab steht.
"""
import sys, os
OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'ps.fonttype': 42})

SUP_F, SUP_E = '#F5DB9E', '#B26A00'   # supervised step
UNS_F, UNS_E = '#9EDECC', '#00795A'   # unsupervised step
NEU_F, NEU_E = '#E9EDF1', '#5A6570'   # data step, no selection
MOD_F, MOD_E = '#9ECAE4', '#0072B2'   # PLS-DA
A1U_F, A1U_E = '#C8EDE0', '#009E73'   # Route A1u
VALA_F, VALA_E = '#F5DB9E', '#B26A00'
VALB_F, VALB_E = '#9EDECC', '#00795A'
INK = '#1a1a1a'

# Reference scale for the text-width measurement. The actual image crop is
# computed at the end from the space actually used (X_CONTENT/H_CONTENT),
# not from these two numbers; they only fix the points-per-unit conversion.
#
# SCALE for print-width compliance: the previous version
# rendered at 292 x 400 mm, well beyond the 180 x 225 mm print area, so the journal
# had to shrink it and the 8-9 pt labels ended up under 5 pt on the page. W_IN and
# every font size below are multiplied by SCALE so the 300 dpi file already IS the
# printed size (225 mm tall, the binding constraint) and no further shrinking happens.
SCALE = 0.5625
W_IN = 8.6 * SCALE
X_MAX = 100.0
PT_PER_UNIT = W_IN * 72.0 / X_MAX
UNIT_IN = W_IN / X_MAX

FS_BODY = 9.6 * SCALE
FS_HEAD = 10.4 * SCALE
FS_LANE = 11.0 * SCALE
FS_PANEL = 15.0 * SCALE
FS_TITLE = 13.0 * SCALE
FS_NOTE = 10.0 * SCALE

PAD_X = 1.6
LINE_H = 3.15
GAP = 2.6


def text_w(s, fs, bold=False):
    prop = FontProperties(family='DejaVu Sans', weight='bold' if bold else 'normal')
    return TextPath((0, 0), s, size=fs, prop=prop).get_extents().width / PT_PER_UNIT


def measure(lines, head):
    w = 0.0
    for i, ln in enumerate(lines):
        w = max(w, text_w(ln, FS_HEAD if (head and i == 0) else FS_BODY, head and i == 0))
    return w + 2 * PAD_X, len(lines) * LINE_H + 2.4


class Canvas:
    def __init__(self):
        self.items, self.texts, self.arrows, self.lines = [], [], [], []

    def box(self, x, y_top, lines, fc, ec, head=True, width=None, height=None):
        w, h = measure(lines, head)
        content_h = h
        if width is not None:
            w = width
        if height is not None:
            h = height
        self.items.append((x, y_top - h, w, h, lines, fc, ec, head, content_h))
        return x, y_top - h, w, h

    def text(self, x, y, s, fs, **kw):
        self.texts.append((x, y, s, fs, kw))

    def arrow(self, x1, y1, x2, y2, color='#5A6570', lw=1.5, ls='-'):
        self.arrows.append((x1, y1, x2, y2, color, lw, ls))

    def line(self, pts, color='#8A94A0', lw=1.1, ls=(0, (4, 3))):
        """Streckenzug ohne Pfeilspitze, fuer Sammel- und Verteilerschienen."""
        self.lines.append((pts, color, lw, ls))

    def draw(self, ax):
        for x, y, w, h, lines, fc, ec, head, content_h in self.items:
            ax.add_patch(FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.10,rounding_size=0.7",
                         linewidth=1.4, edgecolor=ec, facecolor=fc, zorder=2))
            extra = max(0.0, h - content_h)
            cy = y + h - 1.2 - extra / 2.0
            for i, ln in enumerate(lines):
                bold = head and i == 0
                ax.text(x + w / 2, cy, ln, ha='center', va='top',
                        fontsize=FS_HEAD if bold else FS_BODY,
                        fontweight='bold' if bold else 'normal',
                        color=INK, zorder=3)
                cy -= LINE_H
        for pts, color, lw, ls in self.lines:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color,
                    linewidth=lw, linestyle=ls, zorder=1, solid_capstyle='butt')
        for x1, y1, x2, y2, color, lw, ls in self.arrows:
            ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                         mutation_scale=13, linewidth=lw, color=color, zorder=1,
                         shrinkA=0, shrinkB=0, linestyle=ls))
        for x, y, s, fs, kw in self.texts:
            ax.text(x, y, s, fontsize=fs, zorder=3, **kw)


C = Canvas()
LEFT = 3.0


def panel_title(y, label, title):
    C.text(LEFT, y, label, FS_PANEL, fontweight='bold', color=INK, va='center', ha='left')
    C.text(LEFT + text_w(label, FS_PANEL, True) + 1.8, y, title, FS_TITLE,
           fontweight='bold', color=INK, va='center', ha='left')


def chain(y_top, specs, gap=GAP, x0=LEFT + 1.0):
    out, x = [], x0
    for k, (lines, fc, ec, head) in enumerate(specs):
        bx, by, bw, bh = C.box(x, y_top, lines, fc, ec, head)
        if k:
            C.arrow(x - gap, y_top - bh / 2, x - 0.4, y_top - bh / 2)
        out.append((bx, by, bw, bh))
        x = bx + bw + gap
    return out


CONTENT_R = 0.0


def track_right(*edges):
    global CONTENT_R
    CONTENT_R = max(CONTENT_R, *edges)


# ============================ PANEL A ============================
y = -3.0
panel_title(y, '(a)', 'Route A: large data-driven descriptor pool')

y -= 5.4
_dp_lines = ['Descriptor pool', '3,507 descriptors', 'from five libraries']
_cl_lines = ['Cleaning', '2,208 descriptors']
_row_h = max(measure(_dp_lines, True)[1], measure(_cl_lines, True)[1])
_x0 = LEFT + 1.0
_dp = C.box(_x0, y, _dp_lines, NEU_F, NEU_E, True, height=_row_h)
_cl = C.box(_dp[0] + _dp[2] + GAP, y, _cl_lines, NEU_F, NEU_E, True, height=_row_h)
C.arrow(_dp[0] + _dp[2], y - _row_h / 2, _cl[0] - 0.4, y - _row_h / 2)
rowA = [_dp, _cl]
y -= rowA[0][3]

LANES = [
    ('Route A1: unsupervised grouping first (recommended)', UNS_E,
     [(['Correlation grouping of all 2,208', '182 representatives'], UNS_F, UNS_E, False),
      (['Effect size filter on the 182', '35 to 56 pass'], SUP_F, SUP_E, False)]),
    ('Route A1u: label-free choice of group representative', A1U_E,
     [(['Same grouping, representative', 'chosen by centrality'], A1U_F, A1U_E, False),
      (['Effect size filter', '11 to 27 pass'], SUP_F, SUP_E, False)]),
    ('Route A2: supervised filter first (comparison)', SUP_E,
     [(['Effect size filter on all 2,208', '373 to 603 pass'], SUP_F, SUP_E, False),
      (['Correlation grouping', '27 to 50 representatives'], UNS_F, UNS_E, False)]),
]

# Alle drei Spuren bekommen dieselbe Spaltenbreite (Kasten 1 und Kasten 2 je
# ueber alle Routen hinweg), damit sie denselben linken und rechten Rand haben,
# statt je nach Textlaenge unterschiedlich weit zu reichen.
COL_W = [max(measure(specs[ci][0], specs[ci][3])[0] for _, _, specs in LANES) for ci in (0, 1)]

lane_rows = []
A1U_ARROW_Y = None
for li, (label, colour, specs) in enumerate(LANES):
    y -= 4.8
    C.text(LEFT + 1.0, y, label, FS_LANE, style='italic', color=colour, ha='left', va='center')
    y -= 3.0
    row_y = y
    out, x = [], LEFT + 1.0
    for ci, (lines, fc, ec, head) in enumerate(specs):
        bx, by, bw, bh = C.box(x, row_y, lines, fc, ec, head, width=COL_W[ci])
        if ci:
            C.arrow(x - GAP, row_y - bh / 2, x - 0.4, row_y - bh / 2)
            if li == 1:
                A1U_ARROW_Y = row_y - bh / 2
        out.append((bx, by, bw, bh))
        x = bx + bw + GAP
    lane_rows.append(out)
    y -= out[0][3]

# Rechte Spalte: eine gemeinsame Spalte fuer alle PLS-DA-Kaesten in a, b und c,
# so platziert, dass der Pfeil nach unten in jeder Route senkrecht auf die Mitte
# des Kastens darunter (Nested validation required / Simple permutation is
# already honest / VIP-Hinweis) trifft. Diese drei Kaesten liegen bereits alle
# auf derselben Spalte (BUS, Breite NESTED_W), berechnet aus dem Text der
# laengsten der drei, bevor sie gezeichnet werden.
banner = ['Nested validation required',
          'Every label-dependent step, including the choice of a supervised or a label-free',
          'group representative, is repeated inside each leave-one-out fold, and 2,000',
          'permutations scramble the labels through the whole pipeline. An ordinary',
          'permutation with a fixed selection is circular.']
BUS = max(r[-1][0] + r[-1][2] for r in lane_rows) + 5.0
ELBOW_X = BUS - 3.0
NESTED_W = measure(banner, True)[0]
BUS_CENTER = BUS + NESTED_W / 2.0

PLSDA_A = ['PLS-DA', 'one latent variable', 'leave-one-out']
PLSDA_B = ['PLS-DA', 'one latent variable']
PLSDA_C = ['Second PLS-DA', 'reported model']
PLSDA_W = max(measure(PLSDA_A, True)[0], measure(PLSDA_B, True)[0], measure(PLSDA_C, True)[0])
PLSDA_X = BUS_CENTER - PLSDA_W / 2.0
PLSDA_CENTER_X = BUS_CENTER

A1U_TOP = lane_rows[1][1][1] + lane_rows[1][1][3]
# Final set (3 text lines) is taller than the "Effect size filter / 11 to 27 pass"
# box (2 lines) it receives its arrow from, so aligning the two boxes' tops (as
# before) left the incoming arrow above Final set's true centre, and put Final
# set and PLS-DA close enough together that the arrow between them was barely
# visible. Final set is nudged left (more room for that arrow) and both Final
# set and PLS-DA (moved up by the same amount, so the arrow between them stays
# horizontal) are shifted so Final set is centred on the incoming arrow's height.
_finalSet_lines = ['Final set', 'optional top 15', 'by effect size']
_finalSet_w, _finalSet_h = measure(_finalSet_lines, True)
DX_LEFT = 2.0
DY_UP = A1U_ARROW_Y - (A1U_TOP - _finalSet_h / 2.0)
FINAL_TOP = A1U_TOP + DY_UP
finalSet = C.box(PLSDA_X - GAP - _finalSet_w - DX_LEFT, FINAL_TOP, _finalSet_lines, SUP_F, SUP_E, True)
plsdaA = C.box(PLSDA_X, FINAL_TOP, PLSDA_A, MOD_F, MOD_E, True, width=PLSDA_W)
C.arrow(finalSet[0] + finalSet[2], FINAL_TOP - plsdaA[3] / 2, PLSDA_X - 0.4, FINAL_TOP - plsdaA[3] / 2)
rowEnd = [finalSet, plsdaA]

for r in lane_rows:
    xr = r[-1][0] + r[-1][2]
    yr = r[-1][1] + r[-1][3] / 2
    C.line([(xr + 0.5, yr), (ELBOW_X, yr)])
C.line([(ELBOW_X, lane_rows[0][-1][1] + lane_rows[0][-1][3] / 2),
        (ELBOW_X, lane_rows[-1][-1][1] + lane_rows[-1][-1][3] / 2)])
C.arrow(ELBOW_X, A1U_ARROW_Y, finalSet[0] - 0.4, A1U_ARROW_Y, color='#8A94A0', lw=1.1)

# Verteilerschiene links: unter der Bereinigung nach links, dann in die drei Spuren
SPLIT = LEFT - 2.2
_xc = rowA[1][0] + rowA[1][2] / 2
_y_bus = rowA[1][1] - 2.2
C.line([(_xc, rowA[1][1]), (_xc, _y_bus), (SPLIT, _y_bus),
        (SPLIT, lane_rows[-1][0][1] + lane_rows[-1][0][3] / 2)])
for r in lane_rows:
    yr = r[0][1] + r[0][3] / 2
    C.arrow(SPLIT, yr, r[0][0] - 0.4, yr, color='#8A94A0', lw=1.1, ls=(0, (4, 3)))

A2_BOTTOM = lane_rows[2][1][1]
NESTED_Y_TOP = A2_BOTTOM - 3.5
bx, by, bw, bh = C.box(BUS, NESTED_Y_TOP, banner, VALA_F, VALA_E, head=True, width=NESTED_W)

C.arrow(PLSDA_CENTER_X, plsdaA[1], PLSDA_CENTER_X, NESTED_Y_TOP - 0.4, color=VALA_E, lw=1.9)

track_right(rowEnd[1][0] + rowEnd[1][2], bx + bw)
y = by

# ============================ PANEL B ============================
y -= 7.0
panel_title(y, '(b)', 'Route B: small literature-derived descriptor set')
y -= 5.4
rowB_pre = chain(y, [(['Literature set', '25 interpretable', '2D descriptors'], NEU_F, NEU_E, True),
                     (['Correlation grouping only', 'no supervised step at any point',
                       '25 to 12 descriptors'], UNS_F, UNS_E, True)])
plsdaB = C.box(PLSDA_X, y, PLSDA_B, MOD_F, MOD_E, True, width=PLSDA_W)
C.arrow(rowB_pre[-1][0] + rowB_pre[-1][2], y - plsdaB[3] / 2, PLSDA_X - 0.4, y - plsdaB[3] / 2)
rowB = rowB_pre + [plsdaB]
y -= rowB[0][3] + 3.2
bx2, by2, bw2, bh2 = C.box(BUS, y,
                           ['Simple permutation is already honest',
                            'no supervised selection, so no selection step needs nesting'],
                           VALB_F, VALB_E, head=True, width=NESTED_W)
C.arrow(PLSDA_CENTER_X, plsdaB[1], PLSDA_CENTER_X, y - 0.4, color=VALB_E, lw=1.9)
track_right(bx2 + bw2)
y -= bh2

# ============================ PANEL C ============================
y -= 7.0
panel_title(y, '(c)', 'VIP-based selection')
y -= 5.0
C.text(LEFT + 1.0, y, 'The first PLS-DA only ranks descriptors by VIP; the second, fitted on the top 15, is the model whose performance is reported.',
       FS_NOTE, color='#444', ha='left', va='center')
y -= 3.3
y -= 1.0
rowC_pre = chain(y, [(['Route A1 set', '35 to 56 descriptors'], NEU_F, NEU_E, True),
                     (['First PLS-DA', 'VIP ranking only'], MOD_F, MOD_E, True),
                     (['Keep the top 15', 'by VIP'], SUP_F, SUP_E, False)])
plsdaC = C.box(PLSDA_X, y, PLSDA_C, MOD_F, MOD_E, True, width=PLSDA_W)
C.arrow(rowC_pre[-1][0] + rowC_pre[-1][2], y - plsdaC[3] / 2, PLSDA_X - 0.4, y - plsdaC[3] / 2)
rowC = rowC_pre + [plsdaC]
y -= rowC[0][3] + 3.2
bx3, by3, bw3, bh3 = C.box(BUS, y,
                           ['VIP is computed from the labels, so the first model must be nested as well'],
                           VALA_F, VALA_E, head=False, width=NESTED_W)
C.arrow(PLSDA_CENTER_X, plsdaC[1], PLSDA_CENTER_X, y - 0.4, color=VALA_E, lw=1.9)
track_right(bx3 + bw3)
y -= bh3

Y_MIN = y - 15.0
X_CONTENT = CONTENT_R + 3.0

# ============================ ZEICHNEN ============================
H = -Y_MIN
fig = plt.figure(figsize=(X_CONTENT * UNIT_IN, H * UNIT_IN), dpi=300)
# Die Achse muss die ganze Leinwand fuellen, sonst stimmt der Massstab nicht,
# mit dem oben die Kastenbreiten aus der Textbreite berechnet wurden.
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_xlim(-1.0, X_CONTENT)
ax.set_ylim(Y_MIN, 0)
ax.axis('off')
C.draw(ax)

handles = [
    Patch(facecolor=SUP_F, edgecolor=SUP_E, label='Supervised step: uses the activity labels, source of selection leakage'),
    Patch(facecolor=UNS_F, edgecolor=UNS_E, label='Unsupervised step: uses descriptor values only'),
    Patch(facecolor=A1U_F, edgecolor=A1U_E, label='Route A1u: label-free choice of representative'),
    Patch(facecolor=MOD_F, edgecolor=MOD_E, label='PLS-DA model'),
    Patch(facecolor=NEU_F, edgecolor=NEU_E, label='Data step without descriptor selection'),
]
ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.01, 0.0), ncol=2,
          fontsize=FS_NOTE, frameon=False, handlelength=1.6, handleheight=1.2,
          labelspacing=0.55, columnspacing=2.2)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT, 'Fig1_two_route_scheme.' + ext),
                dpi=300, bbox_inches='tight', facecolor='white')
print('Fig. 1 redrawn: %.1f x %.1f in (aspect ratio %.2f), Route B with 12 descriptors.'
      % (X_CONTENT * UNIT_IN, H * UNIT_IN, (H * UNIT_IN) / (X_CONTENT * UNIT_IN)))
