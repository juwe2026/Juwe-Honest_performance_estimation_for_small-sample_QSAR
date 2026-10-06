# -*- coding: utf-8 -*-
"""Rekonstruktion von nested_cv_viptop (Route A2, Auswahl der 15 Deskriptoren
nach VIP statt nach Effektstaerke).

prun_vip.py im abgelegten Archiv ruft ng.nested_cv_viptop auf, aber die Funktion
ist in nested_generic.py nicht vorhanden; sie ist beim Aufraeumen des Codes
verloren gegangen. Die abgelegten Permutationsdateien perm_vip_*.npz existieren
weiter, der Arm ist also gerechnet worden, laesst sich mit dem abgelegten Code
aber nicht nachrechnen.

Hier als nested_cv_a2_viptop neu aufgebaut, Zeile fuer Zeile wie
nested_generic.nested_cv_a2_topk, nur mit _viptop statt der Effektstaerke-
Sortierung in der letzten Auswahlzeile. Die Rekonstruktion wird gegen die
abgelegten Beobachtungswerte geprueft (check_a2vip.py). Der Rueckfallzweig ist
hier bereits korrigiert (0.5 statt ytr.mean()).
"""
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
import nested_generic as ng


def nested_cv_a2_viptop(D, yv, nlv=1, top_k=15, nlv_vip=2):
    N = D['N']; yhat = np.zeros(N)
    for F in D['FOLDS']:
        tr = F['tr']; ytr = yv[tr]; ntr = int(tr.sum())
        n1 = int(ytr.sum())
        if n1 < 1 or n1 == ntr:
            yhat[F['i']] = 0.5; continue
        eff, pv = ng.eff_filter(F['R'], ytr, ntr)
        keep = np.where((eff > 0.3) & (pv <= 0.05))[0]
        if len(keep) < 1:
            yhat[F['i']] = 0.5; continue
        rl = ng.group_reps(F['Xtr'][:, keep], F['R'][:, keep], eff[keep], ntr)
        reps = keep[rl]
        cols = ng._viptop(F['Xtr'], ytr, reps, top_k, nlv_vip) if len(reps) > top_k else reps
        Xtr = F['Xtr'][:, cols]; Xte = F['Xte'][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), ntr - 1), scale=False).fit((Xtr - mu) / sd, ytr)
        yhat[F['i']] = m.predict((Xte - mu) / sd).ravel()[0]
    q2v = 1 - np.sum((yv - yhat) ** 2) / np.sum((yv - yv.mean()) ** 2)
    return q2v, roc_auc_score(yv, yhat), yhat
