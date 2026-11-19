# -*- coding: utf-8 -*-
"""Nested Q² and AUC of the peptides with the descriptor matrix at full precision and rounded to 12, 10 and 8
significant digits (main text, Methods of the second data set, item G10 of 16.11.2026).

Complements input_precision_check.py and rounding_sensitivity.py, which give the number of correlation groups and
the observed hit counts after rounding. Here the nested leave-one-out models themselves are refitted on the rounded
matrices: Routes A1, A1u and A2 with one latent variable, PLS-DA on all 1,275 descriptors and PLS-DA on the 50
label-free representatives without filter. The functions are those of code/ (altorder.alt_nested_cv,
unsup_rep.unsup_nested_cv, nested_generic.nested_cv), the matrix and the activity classes are read as in
revision/lone/run2.py. No permutations are run, so the result gives point estimates only. The values at full
precision reproduce Table 2 of the main text.

Usage:  python rounding_nested_q2.py <archive_root> <empty_work_dir> <out.json>
Nothing is written into the archive. Additional file 3 must be placed in data/ as described in the README, because
code/nested_generic.py reads it on import.
Verified with Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, scikit-learn 1.9.1. Running time about one minute.
"""
import sys, os, json, time
R = sys.argv[1]; os.environ['QSAR_RECALC_DIR'] = sys.argv[2]
os.environ['QSAR_ROUTEA_XLSX'] = R + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'
sys.path.insert(0, os.path.join(R, 'code'))
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
import nested_generic as ng, altorder as ao, unsup_rep as ur

rows = json.load(open(os.path.join(R, 'revision', 'lone', 'rows.json')))
THR = {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}
Y = {k: np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows]) for k in THR}
d0 = pd.read_csv(os.path.join(R, 'revision', 'lone', 'desc_pool.csv'))
DESC = [c for c in d0.columns if c != 'No']; REF = d0[DESC].values.astype(float); N = REF.shape[0]


def sig(A, n):
    A = np.asarray(A, float).copy(); nz = A != 0
    mag = np.zeros_like(A); mag[nz] = np.floor(np.log10(np.abs(A[nz])))
    return np.where(nz, np.round(A / 10.0 ** mag, n - 1) * 10.0 ** mag, 0.0)


def q2(y, yh): return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)


def loo_pls(Xm, y, nlv=1):
    """PLS-DA without selection, autoscaling and the SD < 1e-6 rule inside each training fold."""
    yh = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xt = Xm[tr]; mu = Xt.mean(0); sd = Xt.std(0, ddof=1); ok = sd > 1e-6
        Xt = Xt[:, ok]; mu = mu[ok]; sd = sd[ok]
        m = PLSRegression(n_components=min(nlv, Xt.shape[1], N - 2), scale=False).fit((Xt - mu) / sd, y[tr])
        yh[i] = m.predict((Xm[i:i + 1][:, ok] - mu) / sd).ravel()[0]
    return yh


out = {}
for name, X in [('full', REF), ('sig12', sig(REF, 12)), ('sig10', sig(REF, 10)), ('sig8', sig(REF, 8))]:
    t = time.time(); FOLDS = []
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        FOLDS.append(dict(tr=tr, i=i, Xtr=X[tr], R=rankdata(X[tr], axis=0), Xte=X[i:i + 1]))
    G, _ = ao.build_groups(X); REPS = np.asarray(ur.centrality_reps(X, G))
    out[name] = {'groups': len(G)}
    for k, y in Y.items():
        D = dict(DESC=DESC, X=X, Y=y, N=N, FOLDS=FOLDS)
        arms = [('A1', lambda yy: ao.alt_nested_cv(D, G, yy, 1)[2]), ('A1u', lambda yy: ur.unsup_nested_cv(D, REPS, yy, 1)[2]),
                ('A2', lambda yy: ng.nested_cv(D, yy, 1)[2]), ('all_descriptors', lambda yy: loo_pls(X, yy, 1)),
                ('representatives_no_filter', lambda yy: loo_pls(X[:, REPS], yy, 1))]
        for nm, f in arms:
            yh = f(y); out[name][k + '_' + nm] = dict(q2=float(q2(y, yh)), auc=float(roc_auc_score(y, yh)))
    print(name, 'groups', len(G), '%.0f s' % (time.time() - t), flush=True)
dmax = {}
for lev in ('sig12', 'sig10', 'sig8'):
    for arm in ('A1', 'A1u', 'A2', 'all_descriptors', 'representatives_no_filter'):
        dmax['%s_%s' % (lev, arm)] = max(abs(out[lev][k + '_' + arm]['q2'] - out['full'][k + '_' + arm]['q2']) for k in Y)
json.dump(dict(results=out, largest_abs_change_q2=dmax), open(sys.argv[3], 'w'), indent=1)
print('%-30s %9s %9s %9s %9s' % ('endpoint and model', 'full', 'sig12', 'sig10', 'sig8'))
for k in Y:
    for arm in ('A1', 'A1u', 'A2', 'all_descriptors', 'representatives_no_filter'):
        key = k + '_' + arm
        print('%-30s %9.4f %9.4f %9.4f %9.4f' % (key, *(out[lev][key]['q2'] for lev in ('full', 'sig12', 'sig10', 'sig8'))))
print('largest absolute change of Q2:', json.dumps({k: round(v, 4) for k, v in dmax.items()}))
