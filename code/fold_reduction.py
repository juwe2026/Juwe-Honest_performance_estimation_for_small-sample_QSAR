# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
"""Sensitivitaetsanalyse zu Gutachterpunkt 1.

Die Korrelationsgruppierung von Route B benutzt keine Aktivitaetslabels, sie
benutzt aber die Deskriptorwerte aller 31 Verbindungen. Damit traegt auch die
jeweils ausgelassene Verbindung zur Definition der Merkmale bei. Ob das die
LOO-Schaetzung optimistisch macht, laesst sich nur empirisch klaeren: Hier wird
die Gruppierung in jedem Trainingsfold aus den 30 Trainingsverbindungen neu
bestimmt, einschliesslich der Ausnahmeregel fuer seltene Flags.

Ausgabe: results/recalc/fold_reduction.json (QSAR_RECALC_DIR ueberschreibt das)
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 6 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import sys, json
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import numpy as np, openpyxl, core
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
from scipy.stats import rankdata

THR = 0.7
N = core.N
LIT25 = list(dict.fromkeys(core.LIT25))
assert len(LIT25) == 25, "a priori Satz hat %d Deskriptoren" % len(LIT25)

X25 = core.X_of(LIT25)


def _corr(A):
    """Pearson- und Spearman-Matrix, Spalten mit Nullvarianz bekommen 0."""
    def cc(M):
        M = M - M.mean(0)
        sd = M.std(0)
        ok = sd > 1e-12
        C = np.zeros((M.shape[1], M.shape[1]))
        if ok.sum() > 1:
            sub = np.corrcoef(M[:, ok], rowvar=False)
            idx = np.where(ok)[0]
            C[np.ix_(idx, idx)] = np.nan_to_num(sub)
        return C
    return cc(A), cc(np.column_stack([rankdata(A[:, j]) for j in range(A.shape[1])]))


def _is_binary(v):
    return set(np.unique(v)).issubset({0.0, 1.0})


def group(A, names, thr=THR):
    """Dieselbe Regel wie im Manuskript: seltene Binaerflags (weniger als vier
    Positive) werden vor der Gruppierung ausgenommen, dann greedy nach Grad."""
    P, S = _corr(A)
    L = (np.maximum(np.abs(P), np.abs(S)) > thr)
    np.fill_diagonal(L, False)
    rare = [j for j in range(len(names))
            if _is_binary(A[:, j]) and A[:, j].sum() < 4]
    pool = [j for j in range(len(names)) if j not in rare]
    reps, left = [], list(pool)
    while left:
        deg = {j: int(sum(1 for k in left if k != j and L[j, k])) for j in left}
        best = max(left, key=lambda j: (deg[j], -j))
        partners = [k for k in left if k != best and L[best, k]]
        reps.append(best)
        left = [k for k in left if k != best and k not in partners]
    final = sorted(reps + rare)
    # Die Ausnahme fuer seltene Flags laeuft vor der Gruppierung und laesst exakt
    # redundante Paare durch. Genau dieser Fall tritt global fuer EsterLacton_flag
    # und C_Count auf. Aufloesung wie im Manuskript: das benannte Strukturmerkmal
    # bleibt, der Zaehler faellt weg.
    drop = set()
    for a_i, a in enumerate(final):
        for b in final[a_i + 1:]:
            if a in drop or b in drop:
                continue
            if abs(np.nan_to_num(np.corrcoef(A[:, a], A[:, b])[0, 1])) > 0.99995:
                ab, bb = _is_binary(A[:, a]), _is_binary(A[:, b])
                drop.add(b if (ab and not bb) else a if (bb and not ab) else b)
    return [j for j in final if j not in drop]


# ---------------------------------------------------------------- Gruppierungen
GLOBAL = group(X25, LIT25)
print("Globale Gruppierung: %d Repraesentanten" % len(GLOBAL))
print("  ", ", ".join(LIT25[j] for j in GLOBAL))

FOLD = []
for i in range(N):
    tr = np.ones(N, bool); tr[i] = False
    FOLD.append(group(X25[tr], LIT25))

sizes = [len(f) for f in FOLD]
same = sum(1 for f in FOLD if f == GLOBAL)
print("Foldweise Gruppierung: Groesse %d bis %d, identisch mit global in %d von %d Folds"
      % (min(sizes), max(sizes), same, N))

_union = sorted(set().union(*[set(f) for f in FOLD]))
_inter = sorted(set(FOLD[0]).intersection(*[set(f) for f in FOLD]))
print("  Kern in allen Folds: %d, Vereinigung ueber alle Folds: %d"
      % (len(_inter), len(_union)))
_var = [LIT25[j] for j in _union if j not in _inter]
print("  wechselnde Repraesentanten:", ", ".join(_var) if _var else "keine")


def loo_fold(y, nlv, folds):
    """LOO, in dem die Merkmalsdefinition je Fold aus den Trainingsdaten kommt."""
    yh = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        cols = folds[i]
        Xt = X25[np.ix_(tr, cols)]
        mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd == 0] = 1
        m = PLSRegression(n_components=min(nlv, len(cols), N - 2),
                          scale=False).fit((Xt - mu) / sd, y[tr])
        yh[i] = m.predict((X25[i:i + 1, cols] - mu) / sd).ravel()[0]
    return yh


def run(folds, s, nlv, nperm=2000, seed=702):
    y = core.y_of(s)
    yh = loo_fold(y, nlv, folds)
    oq = core.q2(y, yh); oa = roc_auc_score(y, yh)
    qs = np.empty(nperm); as_ = np.empty(nperm)
    for k in range(nperm):
        yp = np.random.default_rng([seed, k]).permutation(y)
        yhp = loo_fold(yp, nlv, folds)
        qs[k] = core.q2(yp, yhp); as_[k] = roc_auc_score(yp, yhp)
    return dict(q2=float(oq), auc=float(oa),
                p_q2=float((1 + np.sum(qs >= oq - 1e-9)) / (1 + nperm)),
                p_auc=float((1 + np.sum(as_ >= oa - 1e-9)) / (1 + nperm)))


NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
out = {"global_reps": [LIT25[j] for j in GLOBAL],
       "fold_sizes": sizes, "n_identical_folds": same,
       "core": [LIT25[j] for j in _inter], "varying": _var,
       "nperm": NPERM, "results": {}}

for s in core.STRAINS:
    for nlv in (1, 3):
        g = run([GLOBAL] * N, s, nlv, NPERM)
        f = run(FOLD, s, nlv, NPERM)
        out["results"]["%s_%dLV" % (s, nlv)] = {"global": g, "fold": f,
                                                "dq2": f["q2"] - g["q2"],
                                                "dauc": f["auc"] - g["auc"]}
        print("%-7s %dLV  global Q2 %+.3f (p=%.4f)  fold Q2 %+.3f (p=%.4f)  dQ2 %+.3f | "
              "AUC %.3f -> %.3f  dAUC %+.3f"
              % (s, nlv, g["q2"], g["p_q2"], f["q2"], f["p_q2"],
                 f["q2"] - g["q2"], g["auc"], f["auc"], f["auc"] - g["auc"]))

d1 = [v["dq2"] for k, v in out["results"].items() if k.endswith("1LV")]
d3 = [v["dq2"] for k, v in out["results"].items() if k.endswith("3LV")]
out["summary"] = {"mean_dq2_1LV": float(np.mean(d1)), "max_abs_dq2_1LV": float(np.max(np.abs(d1))),
                  "mean_dq2_3LV": float(np.mean(d3)), "max_abs_dq2_3LV": float(np.max(np.abs(d3)))}
print("\nMittlere Differenz Q2 (fold minus global): 1 LV %+.4f, 3 LV %+.4f"
      % (out["summary"]["mean_dq2_1LV"], out["summary"]["mean_dq2_3LV"]))
print("Groesste absolute Differenz: 1 LV %.4f, 3 LV %.4f"
      % (out["summary"]["max_abs_dq2_1LV"], out["summary"]["max_abs_dq2_3LV"]))

json.dump(out, open(_os.path.join(_RES, "fold_reduction.json"), "w"), indent=1)
print("geschrieben: %s" % _os.path.join(_RES, "fold_reduction.json"))
