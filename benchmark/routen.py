"""The three nested routes A1, A1u and A2, vectorised.
Every label-dependent decision is made inside the fold, with the training compounds only.
"""
import numpy as np
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from core_fast import effekt_und_p

P_MAX = 0.05
R_MIN = 0.3


def _pls_vorhersage(Xtr, ytr, xte, sel):
    if len(sel) == 0 or len(set(ytr)) < 2:
        return 0.5
    mu = Xtr[:, sel].mean(0)
    sd = Xtr[:, sel].std(0, ddof=0)
    sd = np.where(sd < 1e-6, 1.0, sd)
    A = (Xtr[:, sel] - mu) / sd
    B = (xte[sel] - mu) / sd
    try:
        pls = PLSRegression(n_components=1, scale=False)
        pls.fit(A, ytr.astype(float))
        return float(np.ravel(pls.predict(B.reshape(1, -1)))[0])
    except Exception:
        return 0.5


def gruppierung(A):
    """Most-connected-first on the matrix A (n x k). Returns groups as index lists."""
    k = A.shape[1]
    if k == 0:
        return []
    with np.errstate(invalid='ignore'):
        P = np.corrcoef(A, rowvar=False)
        S = spearmanr(A).statistic if k > 1 else np.array([[1.0]])
    P = np.atleast_2d(np.nan_to_num(P)); S = np.atleast_2d(np.nan_to_num(S))
    R = (np.abs(P) > 0.7) | (np.abs(S) > 0.7)
    np.fill_diagonal(R, False)
    var = A.var(axis=0, ddof=0)
    offen = set(range(k)); gruppen = []
    while offen:
        idx = np.fromiter(offen, int)
        grade = R[np.ix_(idx, idx)].sum(axis=1)
        mx = grade.max()
        kand = idx[grade == mx]
        hub = int(kand[np.argmax(var[kand])])
        mit = [int(j) for j in idx if j != hub and R[hub, j]]
        gruppen.append([hub] + sorted(mit))
        offen -= set([hub] + mit)
    return gruppen


def fold_a1(Xtr, ytr, xte, gruppen):
    r, p = effekt_und_p(Xtr, ytr)
    reps = [g[int(np.argmax(np.abs(r[g])))] for g in gruppen]
    reps = np.array(reps)
    ok = (p[reps] <= P_MAX) & (np.abs(r[reps]) > R_MIN)
    return _pls_vorhersage(Xtr, ytr, xte, reps[ok])


def fold_a1u(Xtr, ytr, xte, reps_fix):
    r, p = effekt_und_p(Xtr, ytr)
    reps = np.asarray(reps_fix)
    ok = (p[reps] <= P_MAX) & (np.abs(r[reps]) > R_MIN)
    return _pls_vorhersage(Xtr, ytr, xte, reps[ok])


def fold_a2(Xtr, ytr, xte):
    r, p = effekt_und_p(Xtr, ytr)
    durch = np.where((p <= P_MAX) & (np.abs(r) > R_MIN))[0]
    if durch.size == 0:
        return 0.5
    g = gruppierung(Xtr[:, durch])
    reps = np.array([durch[grp[int(np.argmax(np.abs(r[durch[grp]])))]] for grp in g])
    return _pls_vorhersage(Xtr, ytr, xte, reps)


def zentralitaets_reps(A, gruppen):
    """Label-free representatives: largest mean absolute correlation in the group."""
    reps = []
    for g in gruppen:
        if len(g) == 1:
            reps.append(g[0]); continue
        sub = A[:, g]
        with np.errstate(invalid='ignore'):
            C = np.abs(np.nan_to_num(np.corrcoef(sub, rowvar=False)))
        np.fill_diagonal(C, 0.0)
        reps.append(g[int(np.argmax(C.mean(axis=1)))])
    return reps


def loocv(A, y, gruppen, route, reps_fix=None):
    n = len(y)
    pred = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        Xtr, ytr, xte = A[m], y[m], A[i]
        if route == 'A1':
            pred[i] = fold_a1(Xtr, ytr, xte, gruppen)
        elif route == 'A1u':
            pred[i] = fold_a1u(Xtr, ytr, xte, reps_fix)
        else:
            pred[i] = fold_a2(Xtr, ytr, xte)
    ss = ((y - y.mean()) ** 2).sum()
    q2 = 1 - ((y - pred) ** 2).sum() / ss if ss > 0 else 0.0
    return q2, pred
