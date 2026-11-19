# -*- coding: utf-8 -*-
"""Sensitivity of the correlation grouping and of the univariate hit counts to the last digits of the
descriptor values (pseudo-ties).

Descriptor values that are equal in exact arithmetic can differ by floating-point error. The rank-sum
filter and the Spearman correlation operate on ranks, so such pseudo-ties are broken arbitrarily. This
script rounds the descriptor values to a given number of significant digits and reports the number of
correlation groups G and the observed hit counts of Routes A2 and A1
  - for the peptides (revision/lone/desc_pool.csv, read as revision/lone/run2.py reads it),
  - for the monoterpenoids (Additional file 3, placed in data/ as described in the README),
  - and G for the 16 series of the main benchmark analysis (benchmark/datasets_k), with a vectorised
    re-implementation of gruppieren() of benchmark/datasets.py that reproduces all 16 stored G unrounded.

Usage:  python rounding_sensitivity.py <archive_root> <empty_work_dir>
The work directory receives the group caches of code/altorder.py. Nothing is written into the archive.
Verified with Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5. Running time about one minute.
"""
import sys, os, json, glob, numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
R = sys.argv[1]
os.environ['QSAR_RECALC_DIR'] = sys.argv[2]
os.environ['QSAR_ROUTEA_XLSX'] = R + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'
sys.path.insert(0, R + '/code')
import nested_generic as ng, altorder as ao
def sig(A, n):
    A = np.asarray(A, float).copy(); nz = A != 0
    mag = np.zeros_like(A); mag[nz] = np.floor(np.log10(np.abs(A[nz])))
    return np.where(nz, np.round(A / 10.0 ** mag, n - 1) * 10.0 ** mag, 0.0)
def hits(X, G, Y):
    N = X.shape[0]; Rk = rankdata(X, axis=0); out = {}
    for k, y in Y.items():
        e, p = ng.eff_filter(Rk, y, N); ps = (e > 0.3) & (p <= 0.05)
        a1 = np.array([g[int(np.argmax(e[g]))] for g in G]); out[k] = (int(ps.sum()), int(ps[a1].sum()))
    return out
print('=== peptides (archive csv as read by run2.py), rounding level')
rows = json.load(open(R + '/revision/lone/rows.json')); THR = {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}
Y = {k: np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows]) for k in THR}
P = pd.read_csv(R + '/revision/lone/desc_pool.csv'); P = P[[c for c in P.columns if c != 'No']].to_numpy(float)
for n in (None, 15, 14, 13, 12, 10, 8, 6):
    X = P if n is None else sig(P, n); G, _ = ao.build_groups(X); h = hits(X, G, Y)
    print('  digits %-5s G=%d | A2 %s | A1 %s | columns with ties %d' % (n or 'full', len(G), [h[k][0] for k in Y], [h[k][1] for k in Y], int(sum(np.unique(X[:, j]).size < X.shape[0] for j in range(X.shape[1])))))
print('=== monoterpenoids (Additional file 3, 2,208 descriptors)')
D = ng.build('Cleaned_Steps1-4', ycol='S. A. (S. aureus)'); XA = D['X']
YC = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPn': 'S. Pneu (S. pneumoniae)', 'SPy': 'S. pyo (S. pyogenes)'}
YM = {k: (ng._clean[c].astype(str).str.strip() == 'active').astype(float).values for k, c in YC.items()}
for n in (None, 12, 10, 8):
    X = XA if n is None else sig(XA, n); G, _ = ao.build_groups(X); h = hits(X, G, YM)
    print('  digits %-5s G=%d | A2 %s | A1 %s' % (n or 'full', len(G), [h[k][0] for k in YM], [h[k][1] for k in YM]))
print('=== benchmark main analysis: number of groups, stored vs recomputed vs after rounding to 10 digits')
def gruppieren_fast(A):
    Pm = np.corrcoef(A, rowvar=False); S = spearmanr(A).statistic
    Rm = (np.abs(np.nan_to_num(Pm)) > 0.7) | (np.abs(np.nan_to_num(S)) > 0.7); np.fill_diagonal(Rm, False)
    var = A.var(axis=0, ddof=0); offen = np.ones(A.shape[1], bool); ng_ = 0
    while offen.any():
        idx = np.where(offen)[0]; deg = Rm[np.ix_(idx, idx)].sum(1); mx = deg.max(); kand = idx[deg == mx]
        hub = max(kand, key=lambda i: (var[i], -i)); mit = np.where(Rm[hub] & offen)[0]
        offen[hub] = False; offen[mit] = False; ng_ += 1
    return ng_
chg = 0
for f in sorted(glob.glob(R + '/benchmark/datasets_k/DS*_X.csv')):
    ds = os.path.basename(f)[:4]
    if int(ds[2:]) > 16: continue
    A = pd.read_csv(f, float_precision='round_trip').to_numpy(float); g0 = json.load(open(f.replace('_X.csv', '_info.json')))['info']['gruppen']
    g1 = gruppieren_fast(A); g2 = gruppieren_fast(sig(A, 10)); chg += (g2 != g0)
    print('  %s stored G=%3d | recomputed %3d | rounded to 10 digits %3d%s' % (ds, g0, g1, g2, '  <-- differs' if g2 != g0 else ''))
print('  series whose number of groups changes under rounding:', chg, 'of 16')
