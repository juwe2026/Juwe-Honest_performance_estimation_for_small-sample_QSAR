"""Alternative Route-A-Reihenfolge: erst Korrelationsgruppierung (unsupervised),
dann Effektgroessenfilter auf den Repraesentanten."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 4 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import numpy as np, pickle, os
from scipy.stats import rankdata, norm, t as tdist
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
import nested_generic as ng

def _store(X):
    """Cache path keyed on the descriptor matrix itself.

    Without the key the cache returned the groups of whichever matrix was seen
    first, silently, when a second sheet was analysed in the same working
    directory (round 6 report, 9.1).
    """
    import hashlib
    h = hashlib.sha1(np.ascontiguousarray(X).tobytes()).hexdigest()[:8]
    return _os.path.join(_RES, "alt_groups_%s.pkl" % h)

def build_groups(X, thresh=0.7):
    """Unsupervised: Gruppenstruktur aus Deskriptor-Deskriptor-Korrelation (Pearson ODER Spearman).
    Hub = am staerksten vernetzt (Tie-break: hoehere Varianz, rein unsupervised)."""
    n, p = X.shape
    R = rankdata(X, axis=0)
    def adj(M):
        C = np.corrcoef(M, rowvar=False); C = np.nan_to_num(C, nan=0.0); np.fill_diagonal(C, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            tt = C*np.sqrt((n-2)/(1-C**2))
        pv = 2*tdist.sf(np.abs(tt), n-2); pv[np.abs(C) >= 1] = 0.0
        return (np.abs(C) > thresh) & (pv < 0.05)
    ADJ = adj(X) | adj(R); np.fill_diagonal(ADJ, False)
    var = np.nan_to_num(X.var(axis=0))
    remaining = np.ones(p, bool)
    deg = ADJ.sum(1).astype(int)
    groups = []
    while remaining.any():
        d = np.where(remaining, deg, -1)
        mx = d.max()
        if mx <= 0:                      # nur noch Singletons
            for i in np.where(remaining)[0]: groups.append([int(i)])
            break
        cand = np.where(d == mx)[0]
        hub = int(cand[np.argmax(var[cand])])          # unsupervised tie-break
        members = np.concatenate([[hub], np.where(ADJ[hub] & remaining)[0]])
        members = np.unique(members).astype(int)
        groups.append([int(m) for m in members])
        for m in members:                               # Grad inkrementell aktualisieren
            deg -= ADJ[m].astype(int)
        remaining[members] = False
    return groups, ADJ

def get_groups(X):
    STORE = _store(X)
    if os.path.exists(STORE):
        with open(STORE,'rb') as f: return pickle.load(f)
    g, _ = build_groups(X)
    with open(STORE,'wb') as f: pickle.dump(g, f)
    return g

def alt_select(groups, R, ytr, ntr, eff_thr=0.3, p_thr=0.05):
    """Repraesentant je Gruppe = hoechste Effektgroesse (supervised),
    danach Effektgroessenfilter auf den Repraesentanten."""
    eff, pv = ng.eff_filter(R, ytr, ntr)
    reps = np.array([g[int(np.argmax(eff[g]))] for g in groups], dtype=int)
    keep = reps[(eff[reps] > eff_thr) & (pv[reps] <= p_thr)]
    return keep, eff, pv, reps

def alt_nested_cv(D, groups, yv, nlv=1, top_k=None):
    np_ = np; N = D['N']; yhat = np_.zeros(N)
    for F in D['FOLDS']:
        tr = F['tr']; ytr = yv[tr]; ntr = int(tr.sum())
        n1 = int(ytr.sum())
        if n1 < 1 or n1 == ntr: yhat[F['i']] = ytr.mean(); continue
        cols, eff, pv, _ = alt_select(groups, F['R'], ytr, ntr)
        if len(cols) < 1: yhat[F['i']] = ytr.mean(); continue
        if top_k is not None and len(cols) > top_k:
            cols = cols[np_.argsort(-eff[cols])[:top_k]]
        Xtr = F['Xtr'][:, cols]; Xte = F['Xte'][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), ntr-1), scale=False).fit((Xtr-mu)/sd, ytr)
        yhat[F['i']] = m.predict((Xte-mu)/sd).ravel()[0]
    q2 = 1 - np_.sum((yv-yhat)**2)/np_.sum((yv-yv.mean())**2)
    return q2, roc_auc_score(yv, yhat), yhat

def alt_apparent(D, groups, nlv=1, top_k=None):
    """Apparente LOO-CV mit fixer, auf allen Daten bestimmter Auswahl."""
    y = D['Y']; N = D['N']; R = rankdata(D['X'], axis=0)
    cols, eff, pv, reps = alt_select(groups, R, y, N)
    if top_k is not None and len(cols) > top_k:
        cols = cols[np.argsort(-eff[cols])[:top_k]]
    yhat = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xtr = D['X'][tr][:, cols]; Xte = D['X'][i:i+1][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(min(nlv, len(cols), int(tr.sum())-1), scale=False).fit((Xtr-mu)/sd, y[tr])
        yhat[i] = m.predict((Xte-mu)/sd).ravel()[0]
    q2 = 1 - np.sum((y-yhat)**2)/np.sum((y-y.mean())**2)
    return q2, roc_auc_score(y, yhat), len(cols), len(reps)

def alt_viptop(D, groups, yv, nlv=1, top_k=15, nlv_vip=2):
    """A1 + VIP-Auswahl: erstes PLS-DA auf den gefilterten Repraesentanten, Top-k nach VIP, zweites PLS-DA."""
    import numpy as _np
    N=D['N']; yhat=_np.zeros(N)
    for F in D['FOLDS']:
        tr=F['tr']; ytr=yv[tr]; ntr=int(tr.sum()); n1=int(ytr.sum())
        if n1<1 or n1==ntr: yhat[F['i']]=ytr.mean(); continue
        cols,eff,pv,_=alt_select(groups,F['R'],ytr,ntr)
        if len(cols)<1: yhat[F['i']]=ytr.mean(); continue
        if len(cols)>top_k: cols=ng._viptop(F['Xtr'],ytr,cols,top_k,nlv_vip)
        Xtr=F['Xtr'][:,cols]; Xte=F['Xte'][:,cols]
        mu=Xtr.mean(0); sd=Xtr.std(0,ddof=1); sd[sd==0]=1.0
        m=PLSRegression(n_components=min(nlv,len(cols),ntr-1),scale=False).fit((Xtr-mu)/sd,ytr)
        yhat[F['i']]=m.predict((Xte-mu)/sd).ravel()[0]
    q2=1-_np.sum((yv-yhat)**2)/_np.sum((yv-yv.mean())**2)
    return q2, roc_auc_score(yv,yhat), yhat
