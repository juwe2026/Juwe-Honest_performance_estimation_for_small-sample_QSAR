# -*- coding: utf-8 -*-
"""Numerical form of the descriptor pool of the second data set (37 peptides) and its effect on the
correlation grouping and on the observed univariate hit counts.

Many of the peptides are isomers, so many descriptor values are equal in exact arithmetic and differ
only by floating-point error. The rank-sum filter and the Spearman correlation operate on ranks, so
such pseudo-ties are broken by the last digits of the stored values. This script takes the matrix in
several numerical forms and reports, for each, the number of cells that are bit-identical with the
matrix used for all results, the number of columns repeating an earlier column, the number of
correlation groups G and the observed hit counts of Routes A2, A1 and A1u (activity rule A):
  - revision/lone/desc_pool.csv read with the default float parser of pandas, as revision/lone/run2.py
    and run3.py read it (the matrix used for all results),
  - the same file read with float_precision='round_trip' (the decimal values as stored),
  - the sheet Descriptors_1275 of Additional file 6 (optional third argument),
  - the matrix used for all results rounded to 15, 14, 13, 12, 10, 8 and 6 significant digits.

Usage:  python input_precision_check.py <archive_root> <out.json> [<Additional_file_6.xlsx>]
The group caches go to QSAR_RECALC_DIR (set it to an empty directory). Nothing is written into the archive.
Running time below one minute.
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.stats import rankdata
R, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(R, 'code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur

rows = json.load(open(os.path.join(R, 'revision', 'lone', 'rows.json')))
THR = {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}
Y = {k: np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows]) for k in THR}
csv = os.path.join(R, 'revision', 'lone', 'desc_pool.csv')
def mat(df):
    assert list(df.iloc[:, 0].astype(str)) == [r['No'] for r in rows]
    return df.iloc[:, 1:].to_numpy(float)
d0 = pd.read_csv(csv); REF = mat(d0)
forms = [('csv, default parser (matrix used for all results)', REF),
         ('csv, float_precision=round_trip', mat(pd.read_csv(csv, float_precision='round_trip')))]
if len(sys.argv) > 3:
    s = pd.read_excel(sys.argv[3], sheet_name='Descriptors_1275'); assert list(s.columns) == list(d0.columns)
    forms.append(('Additional file 6, sheet Descriptors_1275', mat(s)))
def sig(A, n):
    A = np.asarray(A, float).copy(); nz = A != 0
    mag = np.zeros_like(A); mag[nz] = np.floor(np.log10(np.abs(A[nz])))
    return np.where(nz, np.round(A / 10.0 ** mag, n - 1) * 10.0 ** mag, 0.0)
forms += [('rounded to %d significant digits' % n, sig(REF, n)) for n in (15, 14, 13, 12, 10, 8, 6)]
def repeats(X):
    seen, d = set(), 0
    for j in range(X.shape[1]):
        k = X[:, j].tobytes(); d += k in seen; seen.add(k)
    return d
res = []
for name, X in forms:
    G, _ = ao.build_groups(X); reps = np.asarray(ur.centrality_reps(X, G)); Rk = rankdata(X, axis=0); N = X.shape[0]
    h = {}
    for k, y in Y.items():
        e, p = ng.eff_filter(Rk, y, N); ps = (e > 0.3) & (p <= 0.05)
        a1 = np.array([g[int(np.argmax(e[g]))] for g in G])
        h[k] = [int(ps.sum()), int(ps[a1].sum()), int(ps[reps].sum())]
    ad = np.abs(X - REF); rel = ad / np.maximum(np.maximum(np.abs(X), np.abs(REF)), 1e-300); rel[ad == 0] = 0
    r = dict(form=name, cells=int(X.size), bit_identical_with_reference=int((ad == 0).sum()), max_relative_difference=float(rel.max()),
             columns_repeating_an_earlier_column=repeats(X), G=len(G),
             hits_A2=[h[k][0] for k in Y], hits_A1=[h[k][1] for k in Y], hits_A1u=[h[k][2] for k in Y])
    res.append(r)
    print('%-52s identical %5d of %d  max rel %.2g  repeats %d  G=%d | A2 %s | A1 %s | A1u %s' % (
        name, r['bit_identical_with_reference'], r['cells'], r['max_relative_difference'], r['columns_repeating_an_earlier_column'],
        r['G'], r['hits_A2'], r['hits_A1'], r['hits_A1u']), flush=True)
json.dump(dict(strains=list(Y), versions=dict(numpy=np.__version__, pandas=pd.__version__), forms=res), open(OUT, 'w'), indent=1)
