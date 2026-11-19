"""Run-time measurement for the nested routes. EXCLUSIVELY with permuted labels.
No result with the observed labels is computed or stored.
"""
import json, glob, time, sys
import numpy as np, pandas as pd
from scipy.stats import rankdata, mannwhitneyu, spearmanr
from sklearn.cross_decomposition import PLSRegression

def effektstaerke(x, y):
    n1 = int(y.sum()); n2 = len(y) - n1
    if n1 < 2 or n2 < 2: return 0.0, 1.0
    try: u, p = mannwhitneyu(x[y == 1], x[y == 0], alternative='two-sided')
    except ValueError: return 0.0, 1.0
    return 2 * u / (n1 * n2) - 1, p

def filter_idx(X, y):
    keep = []
    for j in range(X.shape[1]):
        r, p = effektstaerke(X[:, j], y)
        if p <= 0.05 and abs(r) > 0.3: keep.append(j)
    return keep

def gruppen_laden(ds):
    g = json.load(open(f'datasets/{ds}_info.json'))['gruppen']
    return g

def route_a1_fold(Xtr, ytr, Xte, gruppen_idx):
    reps = []
    for g in gruppen_idx:
        best, bs = g[0], -1
        for j in g:
            r, _ = effektstaerke(Xtr[:, j], ytr)
            if abs(r) > bs: bs, best = abs(r), j
        reps.append(best)
    Xr = Xtr[:, reps]
    keep = filter_idx(Xr, ytr)
    if not keep or len(set(ytr)) < 2: return 0.5
    sel = [reps[k] for k in keep]
    mu = Xtr[:, sel].mean(0); sd = Xtr[:, sel].std(0, ddof=0); sd[sd < 1e-6] = 1
    A = (Xtr[:, sel] - mu) / sd; B = (Xte[sel] - mu) / sd
    pls = PLSRegression(n_components=1, scale=False)
    pls.fit(A, ytr.astype(float))
    return float(np.ravel(pls.predict(B.reshape(1, -1)))[0])

def loocv_perm(ds, seed):
    X = pd.read_csv(f'datasets/{ds}_X.csv'); cols = list(X.columns); A = X.values.astype(float)
    y0 = pd.read_csv(f'datasets/{ds}_y.csv')['y'].values
    gi = [[cols.index(c) for c in g if c in cols] for g in gruppen_laden(ds)]
    rng = np.random.default_rng(seed)
    y = rng.permutation(y0)
    pred = np.empty(len(y))
    for i in range(len(y)):
        m = np.ones(len(y), bool); m[i] = False
        pred[i] = route_a1_fold(A[m], y[m], A[i], gi)
    ss = ((y - y.mean()) ** 2).sum()
    return 1 - ((y - pred) ** 2).sum() / ss if ss > 0 else 0.0

if __name__ == '__main__':
    dss = sorted({f.split('/')[-1][:4] for f in glob.glob('datasets/DS*_X.csv')})
    res = {}
    for ds in dss:
        X = pd.read_csv(f'datasets/{ds}_X.csv')
        t0 = time.time(); loocv_perm(ds, 1); dt = time.time() - t0
        res[ds] = dict(n=X.shape[0], p=X.shape[1], sek_pro_permutation=round(dt, 2))
        print(ds, 'n', X.shape[0], 'p', X.shape[1], f'{dt:.1f}s per permutation', flush=True)
    json.dump(res, open('laufzeit.json', 'w'), indent=1)
