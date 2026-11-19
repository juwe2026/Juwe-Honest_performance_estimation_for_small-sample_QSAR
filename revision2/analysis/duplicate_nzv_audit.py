# -*- coding: utf-8 -*-
"""Audit of the archived descriptor matrices: columns whose complete value vector equals that of an earlier column
(exact duplicates, compared bit for bit), and columns that meet the near-zero-variance criterion of Route A
(frequency ratio of the two most common values > 19 and percentage of distinct values < 10). The column sd_le_1e-12
counts columns with a standard deviation of at most 1e-12. In the Route A matrix these are RDF110u and RDF110s, whose
values are of the order 1e-255 for four compounds and zero otherwise, so they are not constant.
Matrices: Route A (Additional file 3, sheet Cleaned_Steps1-4), the peptides (revision/lone/desc_pool.csv) and the
56 benchmark matrices. Usage: python3 duplicate_nzv_audit.py <archive_root> <Additional file 3 .xlsx> <out.json>"""
import sys, glob, os, json
import numpy as np, pandas as pd
A, AF3, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
def audit(X):
    X = X.apply(pd.to_numeric, errors='coerce')
    seen = {}; dup = 0
    for c in X.columns:
        key = X[c].values.astype(float).tobytes()
        if key in seen: dup += 1
        else: seen[key] = c
    nzv = 0; n = len(X)
    for c in X.columns:
        vc = X[c].value_counts(dropna=False).values
        fr = vc[0] / vc[1] if len(vc) > 1 else np.inf
        pu = 100.0 * len(vc) / n
        if fr > 19 and pu < 10: nzv += 1
    small = int((X.std(ddof=0) <= 1e-12).sum())
    return {'n': n, 'p': X.shape[1], 'dup': dup, 'nzv': nzv, 'sd_le_1e-12': small}
out = {}
mt = pd.read_excel(AF3, sheet_name='Cleaned_Steps1-4')
mt = mt.iloc[:, 7:] if mt.shape[1] == 2215 else mt
out['Monoterpenoide'] = audit(mt)
out['Peptide'] = audit(pd.read_csv(os.path.join(A, 'revision/lone/desc_pool.csv'), index_col=0))
for d in ['datasets_k', 'datasets', 'datasets_m']:
    for f in sorted(glob.glob(os.path.join(A, 'benchmark', d, 'DS*_X.csv'))):
        out[d + '/' + os.path.basename(f)[:4]] = audit(pd.read_csv(f))
json.dump(out, open(OUT, 'w'), indent=1)
for k, v in out.items(): print(k, v)
