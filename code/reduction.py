"""Route-B-Korrelationsreduktion auf Basis der PUBLIZIERTEN Matrizen
(Blaetter Pearson_Corr und Spearman_Corr in Zusatzdatei 2)."""
import sys; sys.path.insert(0, '/home/claude/work/recalc')
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

print("Publiziert laut Blatt Reduktion: %d Repraesentanten\n" % len(PUB))
for mode in ("combined", "pearson", "spearman"):
    for ef, lab in ((True, "alt: Ausnahme VOR Gruppierung"), (False, "neu: Gruppierung VOR Ausnahme")):
        reps, cl, ex = reduce(mode, ef)
        mark = ""
        if mode == "combined" and ef:
            mark = "   <-- reproduziert die Publikation" if set(reps) == set(PUB) else "   <-- WEICHT AB"
        print("%-9s | %-32s -> %2d%s" % (mode, lab, len(reps), mark))
print()
reps, cl, ex = reduce("combined", False)
print("Neue Reduktion (combined, korrigierte Reihenfolge): %d Deskriptoren" % len(reps))
for r in reps:
    p = cl.get(r, [])
    print("   %-20s %-16s %s" % (r, "Singleton" if not p else "Cluster (%d)" % (len(p)+1), ", ".join(p)))
