"""SENSITIVITY ANALYSIS MEDIAN RULE (addendum frozen 25.10.2026 21:39, SHA-256 eaf21b7d...): as praed_k.py, on datasets_m/ and ergebnisse_m/.
Originally main analysis (supplementary specification 24.10.2026 17:20, SHA-256 36668775...): as praed.py.
Observed leave-one-out or OOB predictions of all models with the true labels.

Basis for the paired bootstrap (H4) and for sensitivity, specificity, balanced
accuracy and MCC at threshold 0.5. Nothing is decided anew: the same functions as in
run_benchmark.py. As a check, Q2 and AUC must agree exactly with ergebnisse_m/beobachtet.json,
otherwise the script stops.
Usage: python3 praed_m.py   ->  ergebnisse_m/vorhersagen.json
"""
import json
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import run_benchmark_m as rb
import routen

MODELLE = ['app_A1', 'app_A1u', 'app_A2', 'A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter', 'RF']


def vorhersagen(ds, y=None):
    A, y0, g, repsf, fps, _ = rb.lade(ds)
    y = y0 if y is None else y
    out = {}
    out['app_A1'] = rb.pls_loocv(A, y, rb.auswahl_a1(A, y, g))
    out['app_A1u'] = rb.pls_loocv(A, y, rb.auswahl_a1u(A, y, repsf))
    out['app_A2'] = rb.pls_loocv(A, y, rb.auswahl_a2(A, y))
    for r in ('A1', 'A1u', 'A2'):
        out[r] = routen.loocv(A, y, g, r, repsf)[1]
    out['ohneAuswahl'] = rb.pls_loocv(A, y, np.arange(A.shape[1]))
    out['ohneFilter'] = rb.pls_loocv(A, y, np.asarray(repsf))
    rf = RandomForestClassifier(n_estimators=300, min_samples_leaf=1, oob_score=True,
                                bootstrap=True, random_state=0, n_jobs=1).fit(fps, y)
    out['RF'] = np.nan_to_num(rf.oob_decision_function_[:, 1], nan=0.5)
    return y, out


if __name__ == '__main__':
    import glob
    beob = json.load(open('ergebnisse_m/beobachtet.json'))
    dss = sorted({f.split('/')[-1][:4] for f in glob.glob('datasets_m/DS*_X.csv')})
    res = {}
    maxabw = 0.0
    for ds in dss:
        y, P = vorhersagen(ds)
        res[ds] = dict(y=y.tolist(), pred={m: P[m].tolist() for m in MODELLE})
        for m in MODELLE:
            q2, auc = rb.guete(y, P[m])
            s = beob[ds][m]
            abw = max(abs(q2 - s[0]), abs(auc - s[1]))
            maxabw = max(maxabw, abw)
            assert abw < 1e-12, (ds, m, q2, auc, s)
        print(ds, 'ok', flush=True)
    json.dump(res, open('ergebnisse_m/vorhersagen.json', 'w'))
    print('all models reproduced, largest deviation from beobachtet.json:', maxabw)
