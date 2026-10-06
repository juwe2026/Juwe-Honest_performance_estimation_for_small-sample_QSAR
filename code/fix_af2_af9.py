"""Bereinigt zwei Zusatzdateien.

Zusatzdatei 2 enthielt drei Blaetter aus frueheren Rechenlaeufen, deren Zahlen
den Tabellen des Manuskripts widersprechen: Comparison_3_routes mit veralteten
Route-A2-Werten, Table1B_top15 mit p-Werten aus 500 statt 2.000 Permutationen
und einem Namen, der mit Tabelle 1b des Haupttexts kollidiert, sowie
Conclusions mit veraltetem Freitext.

Zusatzdatei 9 trug Blattnamen aus einer frueheren Abbildungsnummerierung
(FigS5/FigS6 statt Fig. S4/Fig. S5).
"""
import json, os, shutil
import numpy as np, openpyxl
from openpyxl.styles import Font

AF2 = "out/Additional_file_2_RouteB_descriptor_data.xlsx"
AF9 = "out/Additional_file_9_source_data_supplementary_figures.xlsx"
SRC9 = "upload/Projekt_Cheminformatics/01_Einreichung/Additional_file_9_source_data_supplementary_figures.xlsx"
AF7 = "out/Additional_file_7_RouteA2_top15_EffSize_vs_VIP.xlsx"

FIN = np.load("repo/QSAR_selection_leakage/results/raw_permutations/final_all.npy",
              allow_pickle=True).item()
RES = json.load(open("recalc/results.json"))
ORDER = ["H.I.", "S.A.", "S.Pneu", "S.Pyo", "P.A."]
FULL = {"H.I.": "H. influenzae", "S.A.": "S. aureus", "S.Pneu": "S. pneumoniae",
        "S.Pyo": "S. pyogenes", "P.A.": "P. aeruginosa"}

# ===================================================================== AF2
wb = openpyxl.load_workbook(AF2)
I = Font(italic=True, size=10)
B = Font(bold=True)

# --- Comparison_3_routes: Route-A2-Spalten aus final_all.npy neu setzen
ws = wb["Comparison_3_routes"]
ws.cell(row=1, column=2, value="(1) Route A2, apparent Q2 (2,208 descriptors)")
ws.cell(row=1, column=3, value="(2) Route A2, nested Q2")
ws.cell(row=1, column=5, value="(3) Route B, Q2 (12 descriptors, nLV = 1, honest)")
for r in range(2, ws.max_row + 1):
    name = ws.cell(row=r, column=1).value
    s = [k for k, v in FULL.items() if v == name]
    if not s:
        continue
    s = s[0]
    f = FIN[s]
    ws.cell(row=r, column=2, value=round(float(f["A2_appq"]), 3))
    ws.cell(row=r, column=3, value=round(float(f["A2_nq"]), 3))
    ws.cell(row=r, column=4, value=round(float(f["A2_p"]), 4))
    b = RES["perm"]["D12"]["1"][s]
    ws.cell(row=r, column=5, value=round(b["q2"], 3))
    ws.cell(row=r, column=6, value=round(b["p_q2"], 4))
n = ws.max_row + 2
for line in [
    "Note. Columns (1) to (3) are Route A2 (supervised filter first, then correlation grouping)",
    "and reproduce the Route A2 rows of main-text Table 1a exactly. Column (5) is Route B with the",
    "12-descriptor modelling set at one latent variable, as in Table 1a. All p-values come from",
    "2,000 whole-pipeline permutations.",
]:
    ws.cell(row=n, column=1, value=line).font = I
    n += 1
print("  Comparison_3_routes an Tabelle 1a angeglichen")

# --- Table1B_top15: p-Werte aus 2.000 Permutationen, Blatt umbenannt
w7 = openpyxl.load_workbook(AF7, data_only=True)["Top15_RouteA2_EffSize_vs_VIP"]
rows7 = {r[0]: r for r in w7.iter_rows(min_row=2, values_only=True) if r[0]}
ws = wb["Table1B_top15"]
for r in range(2, ws.max_row + 1):
    name = ws.cell(row=r, column=1).value
    if name not in rows7:
        continue
    a = rows7[name]           # Strain, appQ2, nestQ2, p(Q2), nestAUC, p(AUC), ...
    ws.cell(row=r, column=3, value=round(float(a[1]), 3))
    ws.cell(row=r, column=5, value=round(float(a[2]), 3))
    ws.cell(row=r, column=7, value=a[3] if isinstance(a[3], str) else round(float(a[3]), 4))
    ws.cell(row=r, column=6, value=round(float(a[4]), 3))
    ws.cell(row=r, column=8, value=a[5] if isinstance(a[5], str) else round(float(a[5]), 4))
    ws.cell(row=r, column=9, value=2000)
n = ws.max_row + 2
for line in [
    "Note. Route A2 restricted to the 15 descriptors with the largest effect size per strain.",
    "These values are identical to Additional file 10 and are NOT the values of main-text Table 3",
    "and Table 4, which report Route A1. Earlier versions of this sheet carried p-values from 500",
    "permutations; they have been replaced by the 2,000-permutation values used throughout.",
]:
    ws.cell(row=n, column=1, value=line).font = I
    n += 1
ws.title = "RouteA2_top15"
print("  Table1B_top15 auf 2.000 Permutationen gesetzt und in RouteA2_top15 umbenannt")

# --- Conclusions: veraltete Zahlen ersetzen
ws = wb["Conclusions"]
FIX = [
    ("S. pyogenes remains weak (Q2 = 0.06, n.s.)",
     "S. pyogenes remains weak (Q2 = 0.05, p = 0.036, nominal only)"),
    ("under Bonferroni (alpha = 0.01)", "under Bonferroni (alpha = 0.0125)"),
    ("H.I. (0.001), S.Pneu (0.0025) and S.A. pass (p(Q2)",
     "H.I. (0.002), S.Pneu (0.001) and S.A. pass (p(Q2)"),
    ("= 0.003; p(AUC) = 0.018 borderline)", "= 0.001; p(AUC) = 0.017 borderline)"),
    ("Q2 = 0.29-0.39, p < 0.01", "Q2 = 0.29-0.39, p not greater than 0.002"),
]
hits = 0
for r in range(1, ws.max_row + 1):
    v = ws.cell(row=r, column=1).value
    if not isinstance(v, str):
        continue
    v2 = v
    for a, b in FIX:
        v2 = v2.replace(a, b)
    if v2 != v:
        ws.cell(row=r, column=1, value=v2)
        hits += 1
n = ws.max_row + 2
for line in [
    "",
    "ADDENDUM (final version). All numbers in this sheet follow main-text Table 1a and Table 1b:",
    "Route B with 12 descriptors at one latent variable gives Q2 = 0.29 / 0.39 / 0.35 / 0.05 / 0.82",
    "with p = 0.002 / 0.001 / 0.001 / 0.036 / <0.001 for H. influenzae, S. aureus, S. pneumoniae,",
    "S. pyogenes and P. aeruginosa. The Bonferroni threshold for the four evaluable endpoints is",
    "alpha = 0.0125. See also sheet Fold_internal_reduction for the sensitivity of these values to",
    "the point at which the unsupervised correlation grouping is performed.",
]:
    ws.cell(row=n, column=1, value=line).font = I
    n += 1
print("  Conclusions: %d Zeilen aktualisiert" % hits)

wb.save(AF2)
print("geschrieben:", AF2)

# ===================================================================== AF9
if not os.path.exists(AF9):
    shutil.copy(SRC9, AF9)
wb9 = openpyxl.load_workbook(AF9)
REN = {"FigS5_LV1_scores": "FigS4_LV1_scores", "FigS6_scores_2comp": "FigS5_scores_2comp"}
for s in list(wb9.sheetnames):
    new = REN.get(s)
    if new is None and s.startswith("S6_"):
        new = "S5_" + s[3:]
    if new and new != s:
        wb9[s].title = new
print("  Zusatzdatei 9, Blattnamen:", wb9.sheetnames)
wb9.save(AF9)
print("geschrieben:", AF9)
