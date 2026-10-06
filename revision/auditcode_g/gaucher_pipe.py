FB=[0]
def _fb():
    FB[0]+=1
    return 0.5

# REVISION OCTOBER 2026: the empty-selection / single-class fallback below previously
# predicted ytr.mean(). Under leave-one-out that equals (S - y_i)/(n - 1) and therefore
# carries the held-out label; it is now the neutral constant 0.5. See CHANGELOG_revision.md.
"""Route A1 / A1u / A2 pipeline applied to the Smit 2007 Gaucher SELDI-TOF data.
Same rules as the monoterpenoid study: correlation grouping (|r| or |rho| > 0.7, p < 0.05),
rank-biserial effect-size filter (|r| > 0.3, p <= 0.05), PLS-DA with ONE latent variable,
autoscaling inside every fold, leave-one-out CV, and whole-pipeline label permutation."""
import numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score, matthews_corrcoef

def load(path):
    df = pd.read_csv(path)
    y = (df['class'].astype(str).str.lower().str.contains('gaucher')).astype(float).values
    X = df.drop(columns=['sample','class']).values.astype(float)
    mz = np.array([float(c) for c in df.drop(columns=['sample','class']).columns])
    return X, y, mz, df['sample'].values

def clean(X, mz):
    keep = np.ones(X.shape[1], bool)
    sd = X.std(0, ddof=1)
    keep &= np.isfinite(sd) & (sd > 1e-12)
    rng = X.max(0) - X.min(0)
    keep &= rng > 1e-12
    return X[:, keep], mz[keep], keep

def eff_filter(R, y, n):
    """rank-biserial effect size + Mann-Whitney p, vectorised over columns."""
    n1 = float(y.sum()); n0 = float(n - n1)
    R1 = R[y == 1].sum(axis=0)
    U1 = R1 - n1*(n1+1)/2.0
    r_rb = 2.0*U1/(n1*n0) - 1.0
    mu = n1*n0/2.0; sig = np.sqrt(n1*n0*(n+1)/12.0)
    Z = (U1 - mu)/sig
    p = 2*norm.sf(np.abs(Z))
    return np.abs(r_rb), p

def build_groups(X, thresh=0.7):
    n, p = X.shape
    R = rankdata(X, axis=0)
    def adj(M):
        C = np.corrcoef(M, rowvar=False); C = np.nan_to_num(C); np.fill_diagonal(C, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            tt = C*np.sqrt((n-2)/(1-C**2))
        pv = 2*tdist.sf(np.abs(tt), n-2); pv[np.abs(C) >= 1] = 0.0
        return (np.abs(C) > thresh) & (pv < 0.05)
    A = adj(X) | adj(R); np.fill_diagonal(A, False)
    var = np.nan_to_num(X.var(axis=0))
    remaining = np.ones(p, bool); deg = A.sum(1).astype(int); groups = []
    while remaining.any():
        d = np.where(remaining, deg, -1); mx = d.max()
        if mx <= 0:
            groups += [[int(i)] for i in np.where(remaining)[0]]; break
        cand = np.where(d == mx)[0]; hub = int(cand[np.argmax(var[cand])])
        mem = np.unique(np.concatenate([[hub], np.where(A[hub] & remaining)[0]])).astype(int)
        groups.append([int(m) for m in mem])
        for m in mem: deg -= A[m].astype(int)
        remaining[mem] = False
    return groups, A

def centrality_reps(X, groups):
    R = rankdata(X, axis=0); reps = []
    for g in groups:
        g = np.asarray(g)
        if len(g) == 1: reps.append(int(g[0])); continue
        Cp = np.corrcoef(X[:, g], rowvar=False); Cs = np.corrcoef(R[:, g], rowvar=False)
        C = np.maximum(np.abs(np.nan_to_num(Cp)), np.abs(np.nan_to_num(Cs)))
        np.fill_diagonal(C, np.nan)
        reps.append(int(g[int(np.nanargmax(np.nanmean(C, axis=1)))]))
    return np.array(sorted(reps), dtype=int)

def group_reps_from(Xk, Rk, effk, n, thresh=0.7):
    """A2: group the survivors, representative = largest effect size."""
    p = Xk.shape[1]
    if p == 0: return np.array([], int)
    def adj(M):
        C = np.corrcoef(M, rowvar=False); C = np.nan_to_num(C); np.fill_diagonal(C, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            tt = C*np.sqrt((n-2)/(1-C**2))
        pv = 2*tdist.sf(np.abs(tt), n-2); pv[np.abs(C) >= 1] = 0.0
        return (np.abs(C) > thresh) & (pv < 0.05)
    A = adj(Xk) | adj(Rk); np.fill_diagonal(A, False)
    remaining = np.ones(p, bool); deg = A.sum(1).astype(int); reps = []
    while remaining.any():
        d = np.where(remaining, deg, -1); mx = d.max()
        if mx <= 0:
            reps += list(np.where(remaining)[0]); break
        cand = np.where(d == mx)[0]; hub = int(cand[np.argmax(effk[cand])])
        mem = np.unique(np.concatenate([[hub], np.where(A[hub] & remaining)[0]])).astype(int)
        reps.append(int(mem[int(np.argmax(effk[mem]))]))
        for m in mem: deg -= A[m].astype(int)
        remaining[mem] = False
    return np.array(sorted(reps), int)

def select(route, X, y, groups, reps_u, eff_thr=0.3, p_thr=0.05):
    n = X.shape[0]; R = rankdata(X, axis=0)
    eff, pv = eff_filter(R, y, n)
    if route == 'A1':
        reps = np.array([g[int(np.argmax(eff[g]))] for g in groups])
        cols = reps[(eff[reps] > eff_thr) & (pv[reps] <= p_thr)]
    elif route == 'A1u':
        cols = reps_u[(eff[reps_u] > eff_thr) & (pv[reps_u] <= p_thr)]
    elif route == 'A2':
        keep = np.where((eff > eff_thr) & (pv <= p_thr))[0]
        if len(keep) == 0: return keep, eff
        rl = group_reps_from(X[:, keep], R[:, keep], eff[keep], n)
        cols = keep[rl]
    return np.sort(cols), eff

def loo_nested(X, y, route, groups, reps_u, nlv=1, top_k=None):
    n = len(y); yh = np.zeros(n)
    for i in range(n):
        tr = np.ones(n, bool); tr[i] = False
        Xtr, ytr = X[tr], y[tr]
        if ytr.sum() < 1 or ytr.sum() == tr.sum(): yh[i] = _fb(); continue
        cols, eff = select(route, Xtr, ytr, groups, reps_u)
        if len(cols) < 1: yh[i] = _fb(); continue
        if top_k is not None and len(cols) > top_k:
            cols = cols[np.argsort(-eff[cols])[:top_k]]
        A = Xtr[:, cols]; B = X[i:i+1, cols]
        mu = A.mean(0); sd = A.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), int(tr.sum())-1), scale=False).fit((A-mu)/sd, ytr)
        yh[i] = m.predict((B-mu)/sd).ravel()[0]
    return yh

def loo_fixed(X, y, cols, nlv=1):
    n = len(y); yh = np.zeros(n)
    for i in range(n):
        tr = np.ones(n, bool); tr[i] = False
        A = X[tr][:, cols]; B = X[i:i+1, cols]
        mu = A.mean(0); sd = A.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), int(tr.sum())-1), scale=False).fit((A-mu)/sd, y[tr])
        yh[i] = m.predict((B-mu)/sd).ravel()[0]
    return yh

def q2(y, yh): return 1 - np.sum((y-yh)**2)/np.sum((y-y.mean())**2)

def metrics(y, yh, thr=0.5):
    p = (yh >= thr).astype(int)
    TP = int(((p==1)&(y==1)).sum()); FN = int(((p==0)&(y==1)).sum())
    FP = int(((p==1)&(y==0)).sum()); TN = int(((p==0)&(y==0)).sum())
    se = TP/(TP+FN) if TP+FN else np.nan; sp = TN/(TN+FP) if TN+FP else np.nan
    try: mcc = matthews_corrcoef(y, p)
    except Exception: mcc = np.nan
    return dict(TP=TP, FN=FN, FP=FP, TN=TN, Sens=se, Spec=sp,
                BA=np.nanmean([se, sp]), MCC=mcc, Acc=(TP+TN)/len(y))
