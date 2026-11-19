"""Switch supplementary file 2 to 12 descriptors. The a priori set of 25 stays unchanged."""
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
import sys, json; sys.path.insert(0,(_ROOT + '/code'))
import openpyxl, core
from openpyxl.styles import Font

R = json.load(open(_WORK + '/results.json'))
CV = json.load(open(_WORK + '/nlvcv.json'))
FULL = {'H.I.':'H. influenzae','S.A.':'S. aureus','S.Pneu':'S. pneumoniae',
        'S.Pyo':'S. pyogenes','P.A.':'P. aeruginosa'}
wb = openpyxl.load_workbook(_WORK + '/AF2_en.xlsx')

def font_of(ws, row=2, colmax=6):
    for c in range(1, colmax+1):
        f = ws.cell(row=row, column=c).font
        if f and f.name: return f.name
    return "Arial"
FONT = font_of(wb["VIP_reduced"])

# ---------------------------------------------------------------- Reduction
ws = wb["Reduction"]
for r in range(2, ws.max_row+1):
    if ws.cell(row=r, column=1).value == "C_Count":
        ws.cell(row=r, column=2, value="Removed: perfectly collinear with EsterLacton_flag")
        ws.cell(row=r, column=3, value="EsterLacton_flag (r = rho = 1.000)")
        ws.cell(row=r, column=4, value=2)
        for c in range(1, 5):
            f = ws.cell(row=r, column=c).font
            ws.cell(row=r, column=c).font = Font(name=f.name or FONT, size=f.size, italic=True)
    if ws.cell(row=r, column=1).value == "EsterLacton_flag":
        ws.cell(row=r, column=2, value="Singleton (rare flag), representative of the C_Count pair")
n = ws.max_row + 2
for line in [
    "Note. The correlation grouping retains 13 representatives. Among them one pair is exactly redundant:",
    "C_Count and EsterLacton_flag have r = rho = 1.000, because Jasmolactone is both the only compound with an",
    "ester or lactone ring and the only one with eleven rather than ten carbon atoms. The pair escaped the grouping",
    "because the exemption for rare binary indicators was applied first. EsterLacton_flag is retained as the",
    "representative, since it names a defined structural feature, whereas C_Count is invariant across the remaining",
    "C10 series. The modelling set therefore comprises 12 descriptors. The a priori set of 25 is unchanged.",
]:
    ws.cell(row=n, column=1, value=line).font = Font(name=FONT, size=10, italic=True)
    n += 1

# ------------------------------------------------------------ PLS-DA_reduced
ws = wb["PLS-DA_reduced"]
ws.cell(row=1, column=3, value="Q2 (nLV = 1, HONEST)")
for r in range(2, ws.max_row+1):
    name = ws.cell(row=r, column=1).value
    s = [k for k, v in FULL.items() if v == name]
    if not s: continue
    s = s[0]
    p1 = R["perm"]["D12"]["1"][s]; cv = CV["D12"][s]
    ws.cell(row=r, column=3, value=round(p1["q2"], 3))
    ws.cell(row=r, column=4, value=round(p1["auc"], 3))
    ws.cell(row=r, column=5, value=round(p1["p_q2"], 4))
    ws.cell(row=r, column=6, value=round(p1["p_auc"], 4))
    ws.cell(row=r, column=7, value=cv["nlv"])
    ws.cell(row=r, column=8, value=round(cv["q2"], 3))
    ws.cell(row=r, column=9, value=round(cv["auc"], 3))
    ws.cell(row=r, column=10, value="%d/%d/%d/%d" % (cv["TP"], cv["FN"], cv["FP"], cv["TN"]))
n = ws.max_row + 2
ws.cell(row=n, column=1, value="12-descriptor modelling set; 2,000 label permutations per model (seed 702, code/run_routeB.py).").font = Font(name=FONT, size=10, italic=True)

# ------------------------------------------------------------- VIP_reduced
ws = wb["VIP_reduced"]
hdr = [ws.cell(row=1, column=c).value for c in range(1, 7)]
keep = [d for d in core.RED12]
for i, d in enumerate(keep):
    ws.cell(row=2+i, column=1, value=d)
    for c in range(2, 7):
        ws.cell(row=2+i, column=c, value=round(R["vip"]["D12"][hdr[c-1]][d], 3))
last = 2 + len(keep)
while ws.max_row >= last:
    ws.delete_rows(last)
ws.cell(row=last+1, column=1, value="C_Count removed (perfectly collinear with EsterLacton_flag); nLV per strain as in VIP_gt1_list.").font = Font(name=FONT, size=10, italic=True)

# -------------------------------------------------------- Comparison_3_routes
ws = wb["Comparison_3_routes"]
ws.cell(row=1, column=5, value="(3) Literature-12 Q2 (nLV = 1, honest)")
for r in range(2, ws.max_row+1):
    name = ws.cell(row=r, column=1).value
    s = [k for k, v in FULL.items() if v == name]
    if not s: continue
    p1 = R["perm"]["D12"]["1"][s[0]]
    ws.cell(row=r, column=5, value=round(p1["q2"], 3))
    ws.cell(row=r, column=6, value=round(p1["p_q2"], 4))

# ---------------------------------------------------------------- Conclusions
ws = wb["Conclusions"]
for r in range(1, ws.max_row+1):
    v = ws.cell(row=r, column=1).value
    if not isinstance(v, str): continue
    v2 = (v.replace("with 13 unsupervised-reduced", "with 12 unsupervised-reduced")
            .replace("13 descriptors combined vs. 14 with Pearson alone",
                     "12 descriptors combined vs. 13 with Pearson alone"))
    if v2 != v: ws.cell(row=r, column=1, value=v2); print("  Conclusions row %d adjusted" % r)
n = ws.max_row + 2
for line in ["", "ADDENDUM: the modelling set comprises 12 descriptors. C_Count was removed because it is",
             "perfectly collinear with EsterLacton_flag (r = rho = 1.000). Q2 changes by at most 0.012,",
             "no significance decision changes, and the confusion matrices are identical. See sheet Reduction."]:
    ws.cell(row=n, column=1, value=line).font = Font(name=FONT, size=10, italic=True); n += 1

wb.save(_WORK + '/Additional_file_2_RouteB_descriptor_data.xlsx')
print("written: out/Additional_file_2_RouteB_descriptor_data.xlsx")
