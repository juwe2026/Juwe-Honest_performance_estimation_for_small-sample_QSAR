# -*- coding: utf-8 -*-
"""How often does the faulty fallback branch actually apply?

The branch (no descriptor passes the filter in the fold, or one class is empty in the
training fold) set the prediction to ytr.mean(). Under LOO this is
(S - y_i)/(n - 1) and therefore reveals the held-out label. What decides
which deposited results are affected is not in which file
the line stands, but whether it is ever executed.

Here, for every arm, we count how many of the 31 folds take the branch with the
true labels, and how many folds and permutations take it under
permuted labels. Counting uses the corrected constant 0.5,
so the counter does not change the result.

Arms in which the branch applies in no permutation are unaffected, nothing
needs to be recomputed there.
"""
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX',
    (_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0, (_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
for _m in ['nested_generic', 'altorder', 'unsup_rep']:
    sys.modules.pop(_m, None)
sys.path.insert(0, (_REV + '/auditcode'))
import nested_generic as ng, altorder as ao, unsup_rep as ur, a2vip, a1uvip
for m in (ng, ao, ur):
    assert (_REV + '/auditcode') in m.__file__, m.__file__
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
p_out = (_OUT + '/fallback_audit.json')
if os.path.exists(p_out):
    out = json.load(open(p_out))
print('permutations per arm:', NP, flush=True)
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
        print('%-26s %-6s true labels %2d/31 folds | permuted %5.2f %% of the folds, '
              '%5.1f %% of the permutations  (%.0fs)'
              % (arm, k, obs_folds, out[key]['folds_permuted_pct'],
                 out[key]['permutations_affected_pct'], time.time() - t), flush=True)
        json.dump(out, open(p_out, 'w'), indent=1)
print('done', flush=True)
