# -*- coding: utf-8 -*-
"""Wie oft greift der fehlerhafte Rueckfallzweig ueberhaupt?

Der Zweig (im Fold passiert kein Deskriptor den Filter, oder eine Klasse ist im
Trainingsfold leer) setzte die Vorhersage auf ytr.mean(). Unter LOO ist das
(S - y_i)/(n - 1) und verraet damit das ausgelassene Label. Entscheidend fuer die
Frage, welche abgelegten Ergebnisse betroffen sind, ist nicht, in welcher Datei
die Zeile steht, sondern ob sie je ausgefuehrt wird.

Hier wird fuer jeden Arm gezaehlt: wie viele der 31 Folds den Zweig bei den
echten Labels nehmen, und wie viele Folds und Permutationen ihn unter
permutierten Labels nehmen. Gezaehlt wird mit der korrigierten Konstanten 0.5,
der Zaehler aendert das Ergebnis also nicht.

Arme, bei denen der Zweig in keiner Permutation greift, sind unberuehrt; dort
braucht nichts nachgerechnet zu werden.
"""
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX',
    '/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0, '/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
for _m in ['nested_generic', 'altorder', 'unsup_rep']:
    sys.modules.pop(_m, None)
sys.path.insert(0, '/home/claude/rev/auditcode')
import nested_generic as ng, altorder as ao, unsup_rep as ur, a2vip, a1uvip
for m in (ng, ao, ur):
    assert '/home/claude/rev/auditcode' in m.__file__, m.__file__
from scipy.stats import rankdata

NP = int(os.environ.get('NP', '200'))
FOLDS = []
for i in range(N):
    tr = np.ones(N, bool); tr[i] = False
    FOLDS.append(dict(tr=tr, i=i, Xtr=X[tr], R=rankdata(X[tr], axis=0), Xte=X[i:i+1]))

ARMS = {
 'A1 (full set)':            (31,  lambda D1, yy: ao.alt_nested_cv(D1, G, yy, 1)[2]),
 'A1u (full set)':           (47,  lambda D1, yy: ur.unsup_nested_cv(D1, REPS_U, yy, 1)[2]),
 'A2 (full set)':            (7,   lambda D1, yy: ng.nested_cv(D1, yy, 1)[2]),
 'A1 top-15 effect size':    (59,  lambda D1, yy: ao.alt_nested_cv(D1, G, yy, 1, 15)[2]),
 'A1 top-15 VIP (1 comp)':   (137, lambda D1, yy: ao.alt_viptop(D1, G, yy, 1, 15, 1)[2]),
 'A2 top-15 effect size':    (11,  lambda D1, yy: ng.nested_cv_a2_topk(D1, yy, 1, 15)[2]),
 'A2 top-15 VIP (2 comp)':   (23,  lambda D1, yy: a2vip.nested_cv_a2_viptop(D1, yy, 1, 15, 2)[2]),
 'A1u top-15 VIP (1 comp)':  (313, lambda D1, yy: a1uvip.viptop_unsup(D1, REPS_U, yy, 1, 15, 1)[2]),
 'A1 top-6 effect size':     (91,  lambda D1, yy: ao.alt_nested_cv(D1, G, yy, 1, 6)[2]),
}
ONLY = os.environ.get('ARMS')
if ONLY:
    ARMS = {k: v for k, v in ARMS.items() if k in ONLY.split('|')}

out = {}
p_out = '/home/claude/rev/out/fallback_audit.json'
if os.path.exists(p_out):
    out = json.load(open(p_out))
print('Permutationen je Arm:', NP, flush=True)
for arm, (seed, fn) in ARMS.items():
    for k in KEYS:
        key = '%s | %s' % (arm, k)
        if key in out:
            continue
        y = Y[k]; D1 = dict(DESC=DESC, X=X, Y=y, N=N, FOLDS=FOLDS)
        t = time.time()
        ng.FB[0] = 0
        fn(D1, y)
        obs_folds = ng.FB[0]
        perm_folds = 0; perm_hit = 0
        for j in range(NP):
            yp = np.random.default_rng([seed, j]).permutation(y)
            ng.FB[0] = 0
            fn(D1, yp)
            perm_folds += ng.FB[0]
            perm_hit += 1 if ng.FB[0] else 0
        out[key] = dict(arm=arm, endpoint=k, seed=seed, nperm=NP,
                        folds_observed=obs_folds,
                        folds_permuted=perm_folds,
                        folds_permuted_pct=round(100.0 * perm_folds / (NP * N), 3),
                        permutations_affected=perm_hit,
                        permutations_affected_pct=round(100.0 * perm_hit / NP, 2),
                        seconds=round(time.time() - t, 1))
        print('%-26s %-6s echte Labels %2d/31 Folds | permutiert %5.2f %% der Folds, '
              '%5.1f %% der Permutationen  (%.0fs)'
              % (arm, k, obs_folds, out[key]['folds_permuted_pct'],
                 out[key]['permutations_affected_pct'], time.time() - t), flush=True)
        json.dump(out, open(p_out, 'w'), indent=1)
print('fertig', flush=True)
