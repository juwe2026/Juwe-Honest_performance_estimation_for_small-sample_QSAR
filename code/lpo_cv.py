# REVISION OCTOBER 2026: the empty-selection / single-class fallback below previously
# predicted ytr.mean(). Under leave-one-out that equals (S - y_i)/(n - 1) and therefore
# carries the held-out label; it is now the neutral constant 0.5. See CHANGELOG_revision.md.
"""Leave-pair-out (LPO) cross-validated c-statistic, as recommended by Geroldinger et al. 2023.
For every (event, non-event) pair the pair is removed, the ENTIRE supervised pipeline is
re-run on the remaining n-2 compounds, and the pair is scored as concordant if the predicted
value of the event exceeds that of the non-event."""
import numpy as np
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
import nested_generic as ng, altorder as ao

def _fit_predict(Xtr, ytr, Xte, cols, nlv=1):
    Xt = Xtr[:, cols]; Xe = Xte[:, cols]
    mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd == 0] = 1.0
    m = PLSRegression(n_components=min(nlv, len(cols), Xt.shape[0]-1), scale=False).fit((Xt-mu)/sd, ytr)
    return m.predict((Xe-mu)/sd).ravel()

def lpo_routeA1(D, groups, nlv=1, top_k=None):
    X = D['X']; y = D['Y']; N = D['N']
    ev = np.where(y == 1)[0]; nev = np.where(y == 0)[0]
    conc = 0.0; total = 0
    for i in ev:
        for j in nev:
            keep = np.ones(N, bool); keep[[i, j]] = False
            Xtr = X[keep]; ytr = y[keep]; ntr = int(keep.sum())
            R = rankdata(Xtr, axis=0)
            cols, eff, pv, _ = ao.alt_select(groups, R, ytr, ntr)
            if len(cols) < 1:
                p = np.array([0.5, 0.5])
            else:
                if top_k is not None and len(cols) > top_k:
                    cols = cols[np.argsort(-eff[cols])[:top_k]]
                p = _fit_predict(Xtr, ytr, X[[i, j]], cols, nlv)
            conc += 1.0 if p[0] > p[1] else (0.5 if p[0] == p[1] else 0.0)
            total += 1
    return conc/total, total

def lpo_fixed(X, y, cols, nlv=1):
    N = len(y); ev = np.where(y == 1)[0]; nev = np.where(y == 0)[0]
    conc = 0.0; total = 0
    for i in ev:
        for j in nev:
            keep = np.ones(N, bool); keep[[i, j]] = False
            p = _fit_predict(X[keep], y[keep], X[[i, j]], cols, nlv)
            conc += 1.0 if p[0] > p[1] else (0.5 if p[0] == p[1] else 0.0)
            total += 1
    return conc/total, total
