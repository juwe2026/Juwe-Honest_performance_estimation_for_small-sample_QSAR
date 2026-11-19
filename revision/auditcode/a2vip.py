# -*- coding: utf-8 -*-
"""Reconstruction of nested_cv_viptop (Route A2, selection of the 15 descriptors
by VIP instead of by effect size).

prun_vip.py in the deposited archive calls ng.nested_cv_viptop, but the function
is not present in nested_generic.py. It was lost when the code was tidied up.
The deposited permutation files perm_vip_*.npz still exist,
so the arm was computed, but it cannot be recomputed with the deposited code.

Rebuilt here as nested_cv_a2_viptop, line by line as
nested_generic.nested_cv_a2_topk, only with _viptop instead of the effect size
sorting in the last selection line. The reconstruction is checked against the
deposited observed values (check_a2vip.py). The fallback branch is
already corrected here (0.5 instead of ytr.mean()).
"""
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
import nested_generic as ng
from nested_generic import _fb


def nested_cv_a2_viptop(D, yv, nlv=1, top_k=15, nlv_vip=2):
    N = D['N']; yhat = np.zeros(N)
    for F in D['FOLDS']:
        tr = F['tr']; ytr = yv[tr]; ntr = int(tr.sum())
        n1 = int(ytr.sum())
        if n1 < 1 or n1 == ntr:
            yhat[F['i']] = _fb(); continue
        eff, pv = ng.eff_filter(F['R'], ytr, ntr)
        keep = np.where((eff > 0.3) & (pv <= 0.05))[0]
        if len(keep) < 1:
            yhat[F['i']] = _fb(); continue
        rl = ng.group_reps(F['Xtr'][:, keep], F['R'][:, keep], eff[keep], ntr)
        reps = keep[rl]
        cols = ng._viptop(F['Xtr'], ytr, reps, top_k, nlv_vip) if len(reps) > top_k else reps
        Xtr = F['Xtr'][:, cols]; Xte = F['Xte'][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), ntr - 1), scale=False).fit((Xtr - mu) / sd, ytr)
        yhat[F['i']] = m.predict((Xte - mu) / sd).ravel()[0]
    q2v = 1 - np.sum((yv - yhat) ** 2) / np.sum((yv - yv.mean()) ** 2)
    return q2v, roc_auc_score(yv, yhat), yhat
