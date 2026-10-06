"""Route A mit VOLLSTAENDIG unsupervidierter Korrelationsgruppierung:
Repraesentant = Mitglied mit hoechster Zentralitaet (mittlerer |r| innerhalb der Gruppe),
also OHNE jede Verwendung der Labels. Danach erst der Effektgroessenfilter."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 6 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import numpy as np, pickle, os
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
import nested_generic as ng, altorder as ao

def _repstore(X):
    """Cache path keyed on the descriptor matrix itself; see altorder._store."""
    import hashlib
    h = hashlib.sha1(np.ascontiguousarray(X).tobytes()).hexdigest()[:8]
    return _os.path.join(_RES, "unsup_reps_%s.pkl" % h)

def centrality_reps(X, groups):
    """Repraesentant je Gruppe = hoechste mittlere |Korrelation| zu den anderen Mitgliedern.
    Nutzt ausschliesslich Deskriptorwerte, keine Labels."""
    REPSTORE = _repstore(X)
    if os.path.exists(REPSTORE):
        with open(REPSTORE,'rb') as f: return pickle.load(f)
    R = rankdata(X, axis=0)
    reps = []
    for g in groups:
        g = np.asarray(g)
        if len(g) == 1:
            reps.append(int(g[0])); continue
        Cp = np.corrcoef(X[:, g], rowvar=False); Cs = np.corrcoef(R[:, g], rowvar=False)
        C = np.maximum(np.abs(np.nan_to_num(Cp)), np.abs(np.nan_to_num(Cs)))
        np.fill_diagonal(C, np.nan)
        cent = np.nanmean(C, axis=1)
        reps.append(int(g[int(np.nanargmax(cent))]))
    reps = np.array(sorted(reps), dtype=int)
    with open(REPSTORE,'wb') as f: pickle.dump(reps, f)
    return reps

def unsup_select(reps, R, ytr, ntr, eff_thr=0.3, p_thr=0.05):
    eff, pv = ng.eff_filter(R, ytr, ntr)
    keep = reps[(eff[reps] > eff_thr) & (pv[reps] <= p_thr)]
    return keep, eff, pv

def unsup_nested_cv(D, reps, yv, nlv=1, top_k=None):
    N = D['N']; yhat = np.zeros(N)
    for F in D['FOLDS']:
        tr = F['tr']; ytr = yv[tr]; ntr = int(tr.sum())
        n1 = int(ytr.sum())
        if n1 < 1 or n1 == ntr: yhat[F['i']] = 0.5; continue
        cols, eff, pv = unsup_select(reps, F['R'], ytr, ntr)
        if len(cols) < 1: yhat[F['i']] = 0.5; continue
        if top_k is not None and len(cols) > top_k:
            cols = cols[np.argsort(-eff[cols])[:top_k]]
        Xtr = F['Xtr'][:, cols]; Xte = F['Xte'][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), ntr-1), scale=False).fit((Xtr-mu)/sd, ytr)
        yhat[F['i']] = m.predict((Xte-mu)/sd).ravel()[0]
    q2 = 1 - np.sum((yv-yhat)**2)/np.sum((yv-yv.mean())**2)
    return q2, roc_auc_score(yv, yhat), yhat
