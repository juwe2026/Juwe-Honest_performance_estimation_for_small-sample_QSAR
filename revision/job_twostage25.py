# -*- coding: utf-8 -*-
"""G2-20: Table S8 (zweistufiger VIP>1-Ansatz auf den 25 Literaturdeskriptoren) und VIP_gt1_list mit NumStereoCenters
von 1,8-Cineol = 0. Algorithmus unveraendert aus code/lit_plsda.py (nLV = Q2-optimale Zahl des 12er-Satzes, wie Table S8).
Zuerst Reproduktion der publizierten Werte mit Cineol = 2."""
import sys, json, numpy as np
sys.path.insert(0, '/home/claude/rev/QSAR_selection_leakage/code')
import core
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
def asc(Xtr, Xte):
    mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1; return (Xtr - mu) / sd, (Xte - mu) / sd
def loo(X, y, nlv, selector=None):
    n = len(y); yh = np.zeros(n)
    for i in range(n):
        tr = np.ones(n, bool); tr[i] = False; Xtr, Xte = X[tr], X[i:i+1]; ytr = y[tr]
        cols = selector(Xtr, ytr, nlv) if selector else np.arange(X.shape[1])
        if len(cols) < 1: yh[i] = ytr.mean(); continue
        a, b = asc(Xtr[:, cols], Xte[:, cols]); use = min(nlv, len(cols), tr.sum() - 1)
        yh[i] = PLSRegression(use, scale=False).fit(a, ytr).predict(b).ravel()[0]
    return yh
def q2(y, yh): return 1 - np.sum((y - yh)**2) / np.sum((y - y.mean())**2)
def vip(m, X):
    T = m.x_scores_; W = m.x_weights_; Q = m.y_loadings_.ravel(); A = T.shape[1]
    ssy = np.array([(Q[a]**2) * (T[:, a]**2).sum() for a in range(A)]); p = X.shape[1]
    return np.array([np.sqrt(p * np.sum(ssy * (W[j, :] / np.linalg.norm(W, axis=0))**2) / ssy.sum()) for j in range(p)])
def fit_vip(X, y, nlv):
    mu = X.mean(0); sd = X.std(0, ddof=1); sd[sd == 0] = 1; Xs = (X - mu) / sd
    return vip(PLSRegression(nlv, scale=False).fit(Xs, y), Xs)
def sel(Xtr, ytr, nlv):
    V = fit_vip(Xtr, ytr, min(nlv, Xtr.shape[1])); c = np.where(V > 1)[0]; return c if len(c) else np.arange(Xtr.shape[1])
L25 = sorted(set(core.LIT25)); X2 = core.X_of(L25); ci = core.COMPOUNDS.index('1,8-Cineole'); j = L25.index('NumStereoCenters')
X0 = X2.copy(); X0[ci, j] = 0
NLV = {'H.I.': 4, 'S.A.': 4, 'S.Pneu': 3, 'S.Pyo': 1, 'P.A.': 4}
PUB = {'H.I.': (0.287, 0.810, 0.161, 0.805, 8), 'S.A.': (0.498, 0.931, 0.424, 0.877, 8), 'S.Pneu': (0.525, 0.960, 0.305, 0.853, 8),
       'S.Pyo': (0.318, 0.868, 0.245, 0.841, 10), 'P.A.': (0.998, 1.000, 0.995, 1.000, 6)}
out = {}
for s, nlv in NLV.items():
    y = core.y_of(s); r = {}
    for tag, X in (('cineole2', X2), ('cineole0', X0)):
        V = fit_vip(X, y, nlv); nc = np.where(V > 1)[0]
        yn = loo(X[:, nc], y, nlv); yg = loo(X, y, nlv, selector=sel)
        r[tag] = dict(q2_naive=q2(y, yn), auc_naive=roc_auc_score(y, yn), q2_nested=q2(y, yg), auc_nested=roc_auc_score(y, yg),
                      n_vip_gt1=int(len(nc)), vip_gt1=[L25[k] for k in nc[np.argsort(-V[nc])]], vip=dict(zip(L25, V.tolist())))
    got = tuple(round(r['cineole2'][k], 3) for k in ('q2_naive', 'auc_naive', 'q2_nested', 'auc_nested')) + (r['cineole2']['n_vip_gt1'],)
    print(s, nlv, 'reproduced' if got == PUB[s] else 'MISMATCH %s vs %s' % (got, PUB[s]),
          '| cineole0', tuple(round(r['cineole0'][k], 3) for k in ('q2_naive', 'auc_naive', 'q2_nested', 'auc_nested')), r['cineole0']['n_vip_gt1'],
          '| VIP>1 lost/gained', set(r['cineole2']['vip_gt1']) - set(r['cineole0']['vip_gt1']), set(r['cineole0']['vip_gt1']) - set(r['cineole2']['vip_gt1']), flush=True)
    out[s] = dict(nlv=nlv, **r)
json.dump(out, open('/home/claude/rev/out/twostage25_cineole0.json', 'w'), indent=1)
