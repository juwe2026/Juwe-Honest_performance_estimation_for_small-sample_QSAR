# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
"""Ordinary (non-nested, 'circular') permutation test: descriptors are selected ONCE on the
true labels (the 'apparent' set), then only the PLS-DA model is refit under permuted labels
via LOO-CV using that FIXED descriptor set. This is the naive test that the paper argues is
optimistically biased; it is computed here deliberately, for contrast with the honest nested test."""
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

def loo_fixed(X, y, cols, nlv=1):
    n = len(y); yh = np.zeros(n)
    for i in range(n):
        tr = np.ones(n, bool); tr[i] = False
        Xt = X[tr][:, cols]; Xe = X[i:i+1][:, cols]
        mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd == 0] = 1
        m = PLSRegression(n_components=min(nlv, len(cols), int(tr.sum())-1), scale=False).fit((Xt-mu)/sd, y[tr])
        yh[i] = m.predict((Xe-mu)/sd).ravel()[0]
    return yh

def q2(y, yh):
    return 1 - np.sum((y-yh)**2) / np.sum((y-y.mean())**2)

def ordinary_perm_test(X, y, cols, nlv=1, n_perm=1000, seed=101):
    """Returns obs_q2, obs_auc, p_q2, p_auc using a FIXED descriptor set across all permutations."""
    yh_obs = loo_fixed(X, y, cols, nlv)
    oq = q2(y, yh_obs); oa = roc_auc_score(y, yh_obs)
    rng = np.random.default_rng(seed)
    q2s = np.zeros(n_perm); aucs = np.zeros(n_perm)
    for k in range(n_perm):
        yp = rng.permutation(y)
        yh = loo_fixed(X, yp, cols, nlv)
        q2s[k] = q2(yp, yh); aucs[k] = roc_auc_score(yp, yh)
    pq = (1 + np.sum(q2s >= oq - 1e-9)) / (1 + n_perm)
    pa = (1 + np.sum(aucs >= oa - 1e-9)) / (1 + n_perm)
    return oq, oa, pq, pa, q2s, aucs
