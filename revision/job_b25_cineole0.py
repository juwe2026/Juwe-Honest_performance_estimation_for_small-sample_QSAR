# -*- coding: utf-8 -*-
"""G2-20: Route B with all 25 descriptors, NumStereoCenters of 1,8-cineole = 0 (CIP-based).
Same permutation scheme as the archive permB_25_*.npz (code/run_routeB.py: default_rng([701, k]), 2,000 permutations,
LOO, PLS-DA with 1 LV). First it is checked that the scheme with the old value (cineole = 2) reproduces the archived null values,
then it is recomputed with cineole = 0. p with treatment of ties (pval_tie.py)."""
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
import sys, json, numpy as np
from multiprocessing import Pool
from common import *
import core
from sklearn.metrics import roc_auc_score
RP = (_ROOT + '/results/raw_permutations/')
L25 = sorted(set(core.LIT25)); assert len(L25) == 25
X25 = core.X_of(L25); ci = core.COMPOUNDS.index('1,8-Cineole'); j = L25.index('NumStereoCenters')
assert X25[ci, j] == 2
X0 = X25.copy(); X0[ci, j] = 0
S = {'HI': 'H.I.', 'SA': 'S.A.', 'SPneu': 'S.Pneu', 'SPyo': 'S.Pyo', 'PA': 'P.A.'}
def one(args):
    k, X, kk = args; y = core.y_of(S[k]); yp = np.random.default_rng([701, kk]).permutation(y); h = loo_pls(X, yp, 1)
    return q2(yp, h), roc_auc_score(yp, h)
if __name__ == '__main__':
    out = {}
    with Pool(8) as P:
        for k in S:
            z = np.load(RP + 'permB_25_%s.npz' % k); y = core.y_of(S[k])
            chk = P.map(one, [(k, X25, kk) for kk in range(10)])
            dq = max(abs(a - b) for (a, _), b in zip(chk, z['q2s'][:10])); assert dq < 1e-9, (k, dq)
            assert abs(q2(y, loo_pls(X25, y, 1)) - float(z['obs_q2'])) < 1e-9
            r = P.map(one, [(k, X0, kk) for kk in range(2000)]); qs = np.array([a for a, _ in r]); as_ = np.array([b for _, b in r])
            h = loo_pls(X0, y, 1); oq, oa = q2(y, h), roc_auc_score(y, h)
            out[k] = dict(q2=oq, auc=oa, p_q2=pval(qs, oq), p_auc=pval(as_, oa), n_perm=2000, seed=701, **mcc_conf(y, h),
                          q2_cineole2=float(z['obs_q2']), auc_cineole2=float(z['obs_auc']),
                          p_q2_cineole2=pval(z['q2s'], float(z['obs_q2'])), p_auc_cineole2=pval(z['aucs'], float(z['obs_auc'])),
                          scheme_check_max_dq=float(dq))
            np.savez_compressed((_OUT + '/permB_25_cineole0_%s.npz') % k, q2s=qs, aucs=as_, obs_q2=oq, obs_auc=oa)
            print(k, {a: round(b, 4) for a, b in out[k].items() if isinstance(b, float)}, flush=True)
    json.dump(out, open((_OUT + '/B25_cineole0_seed701.json'), 'w'), indent=1)
