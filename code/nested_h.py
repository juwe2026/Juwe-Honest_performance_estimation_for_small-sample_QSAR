import numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

# --- Daten laden (bereinigte 2208er-Basis) ---
big = pd.ExcelFile('/mnt/user-data/outputs/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
clean = pd.read_excel(big, 'Cleaned_Steps1-4')
hi = pd.read_excel(big, 'CS1-4_H.I.')
DESC = [c for c in clean.columns if '::' in c and c in hi.columns]
X = hi[DESC].apply(pd.to_numeric, errors='coerce').values.astype(float)
Y = (hi['H. I. (H. influenzae)'].astype(str).str.strip() == 'active').astype(float).values
N = len(Y)

# --- Precompute: pro LOO-Fold die Trainings-Rangmatrix (haengt nur an X, nicht an y) ---
FOLDS = []
for i in range(N):
    tr = np.ones(N, bool); tr[i] = False
    Xtr = X[tr]
    R = rankdata(Xtr, axis=0)            # Spearman-Basis + MW-Raenge
    FOLDS.append(dict(tr=tr, i=i, Xtr=Xtr, R=R, Xte=X[i:i+1]))

def eff_filter(R, ytr, ntr):
    n1 = int(ytr.sum()); n0 = ntr - n1
    U1 = R[ytr.astype(bool)].sum(0) - n1*(n1+1)/2.0
    rrb = 2*U1/(n1*n0) - 1.0
    mu = n1*n0/2.0; sigma = np.sqrt(n1*n0*(ntr+1)/12.0)
    z = (np.abs(U1-mu)-0.5)/sigma
    p = 2*norm.sf(z)
    return np.abs(rrb), p

def group_reps(Xk, Rk, effk, ntr):
    k0 = Xk.shape[1]
    if k0 <= 1:
        return np.arange(k0, dtype=int)
    # Korrelationen (Pearson auf Werten, Spearman auf Raengen) + t-Approx-p
    def cor_p(M):
        C = np.corrcoef(M, rowvar=False); C = np.nan_to_num(C, nan=0.0)
        np.fill_diagonal(C, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            tt = C*np.sqrt((ntr-2)/(1-C**2))
        p = 2*tdist.sf(np.abs(tt), ntr-2); p[np.abs(C) >= 1] = 0.0
        return C, p
    P, pP = cor_p(Xk); S, pS = cor_p(Rk)
    ADJ = ((np.abs(P) > 0.7) & (pP < 0.05)) | ((np.abs(S) > 0.7) & (pS < 0.05))
    np.fill_diagonal(ADJ, False)
    k = Xk.shape[1]; remaining = set(range(k)); reps = []
    partners = {j: set(np.where(ADJ[j])[0]) for j in range(k)}
    while remaining:
        cnt = {j: len(partners[j] & remaining) for j in remaining}
        mx = max(cnt.values()) if cnt else 0
        if mx == 0:
            reps.extend(sorted(remaining)); break
        hub = sorted([j for j in remaining if cnt[j] == mx], key=lambda j: (-effk[j], j))[0]
        members = [hub] + sorted(partners[hub] & remaining, key=lambda j: (-effk[j], j))
        rep = max(members, key=lambda j: effk[j])       # Repraesentant = hoechste Effektgroesse
        reps.append(rep); remaining -= set(members)
    return np.array(sorted(reps), dtype=int)

def nested_cv(yv, nlv=1):
    yhat = np.zeros(N); nfeat = []
    for F in FOLDS:
        tr = F['tr']; ytr = yv[tr]; ntr = int(tr.sum())
        n1 = int(ytr.sum()); n0 = ntr - n1
        if n1 < 1 or n0 < 1:
            yhat[F['i']] = ytr.mean(); nfeat.append(0); continue
        eff, p = eff_filter(F['R'], ytr, ntr)
        keep = np.where((eff > 0.3) & (p <= 0.05))[0]
        if len(keep) < 1:
            yhat[F['i']] = ytr.mean(); nfeat.append(0); continue
        reps_local = group_reps(F['Xtr'][:, keep], F['R'][:, keep], eff[keep], ntr)
        cols = keep[reps_local]; nfeat.append(len(cols))
        Xtr = F['Xtr'][:, cols]; Xte = F['Xte'][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        use = min(nlv, len(cols), ntr-1)
        m = PLSRegression(n_components=use, scale=False).fit((Xtr-mu)/sd, ytr)
        yhat[F['i']] = m.predict((Xte-mu)/sd).ravel()[0]
    q2 = 1 - np.sum((yv-yhat)**2)/np.sum((yv-yv.mean())**2)
    auc = roc_auc_score(yv, yhat)
    return q2, auc, yhat, nfeat
