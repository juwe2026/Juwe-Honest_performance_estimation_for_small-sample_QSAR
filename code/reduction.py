"""Route B correlation reduction on the basis of the PUBLISHED matrices
(sheets Pearson_Corr and Spearman_Corr in supplementary file 2)."""
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
import sys; sys.path.insert(0, (_ROOT + '/code'))
import numpy as np, openpyxl, core

THR = 0.7
wb = openpyxl.load_workbook(core.SRC, data_only=True)

def matrix(sheet):
    rows = list(wb[sheet].values); hdr = list(rows[0])[1:]
    names = [r[0] for r in rows[1:] if r[0]]
    M = {}
    for r in rows[1:]:
        if not r[0]: continue
        for k, n in enumerate(hdr):
            v = r[k+1]
            if isinstance(v, (int, float)): M[(r[0], n)] = float(v)
    return names, M

NAMES, P = matrix("Pearson_Corr")
_, S = matrix("Spearman_Corr")

PUB = {r[0]: (r[1], [x.strip() for x in str(r[2]).split(",")] if r[2] else [])
       for r in list(wb["Reduktion"].values)[1:] if r[0]}
RARE = [n for n, (typ, _) in PUB.items() if "seltenes Flag" in str(typ)]

def adj(mode):
    def link(a, b):
        p = abs(P.get((a, b), P.get((b, a), 0.0))); s = abs(S.get((a, b), S.get((b, a), 0.0)))
        return {"combined": max(p, s), "pearson": p, "spearman": s}[mode] > THR
    return link

def reduce(mode, exempt_first):
    link = adj(mode)
    pool = [n for n in NAMES if not (exempt_first and n in RARE)]
    exempt = [n for n in NAMES if exempt_first and n in RARE]
    reps, clusters, left = [], {}, list(pool)
    while left:
        deg = {n: sum(1 for m in left if m != n and link(n, m)) for n in left}
        best = max(left, key=lambda n: (deg[n], -NAMES.index(n)))
        partners = [m for m in left if m != best and link(best, m)]
        reps.append(best); clusters[best] = partners
        left = [m for m in left if m != best and m not in partners]
    return reps + exempt, clusters, exempt

print("Published according to sheet Reduction: %d representatives\n" % len(PUB))
for mode in ("combined", "pearson", "spearman"):
    for ef, lab in ((True, "old: exemption BEFORE grouping"), (False, "new: grouping BEFORE exemption")):
        reps, cl, ex = reduce(mode, ef)
        mark = ""
        if mode == "combined" and ef:
            mark = "   <-- reproduces the publication" if set(reps) == set(PUB) else "   <-- DIFFERS"
        print("%-9s | %-32s -> %2d%s" % (mode, lab, len(reps), mark))
print()
reps, cl, ex = reduce("combined", False)
print("New reduction (combined, corrected order): %d descriptors" % len(reps))
for r in reps:
    p = cl.get(r, [])
    print("   %-20s %-16s %s" % (r, "Singleton" if not p else "Cluster (%d)" % (len(p)+1), ", ".join(p)))
