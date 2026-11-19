# Revision (EW 03.10.2026): Table S5 (continued) of Revision 1 for Route A1 and Route B, 1-4 LV.
# R2 = fit to all 31 compounds (Route A1: selection on all compounds), as Table S10 of Revision 1,
# Q2 apparent = selection once on all compounds, then LOO (Route A1, as Table 1),
# Q2 nested = selection repeated in every fold (Route A1), Route B: fixed 12 descriptors, LOO.
# Deterministic, no permutations. Same code as jobI2_fallback.py (fixcode).
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
import sys, json
sys.path.insert(0, (_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
for _m in ['nested_generic', 'altorder', 'unsup_rep', 'core']: sys.modules.pop(_m, None)
sys.path.insert(0, (_ROOT + '/code'))
import altorder as ao
assert (_ROOT + '/code') in ao.__file__
sys.path.insert(1, (_ROOT + '/code'))
import core
FOLDS = []
for i in range(N):
    tr = np.ones(N, bool); tr[i] = False
    FOLDS.append(dict(tr=tr, i=i, Xtr=X[tr], R=rankdata(X[tr], axis=0), Xte=X[i:i + 1]))
def fit_all(Xm, y, L):
    mu = Xm.mean(0); sd = Xm.std(0, ddof=1); sd[sd == 0] = 1.0; Xs = (Xm - mu) / sd
    m = PLSRegression(min(L, Xm.shape[1], len(y) - 1), scale=False).fit(Xs, y)
    return m.predict(Xs).ravel()
S = {'HI': 'H.I.', 'SA': 'S.A.', 'SPneu': 'S.Pneu', 'SPyo': 'S.Pyo', 'PA': 'P.A.'}
X12 = core.X_of(core.RED12)
out = {'A1': {}, 'B': {}}
for k in KEYS:
    y = Y[k]; D1 = dict(DESC=DESC, X=X, Y=y, N=N, FOLDS=FOLDS)
    cols, eff, pv, reps = ao.alt_select(G, rankdata(X, axis=0), y, N)
    nsel_fold = [len(ao.alt_select(G, F['R'], y[F['tr']], N - 1)[0]) for F in FOLDS]
    rA = {'n_selected_all': int(len(cols)), 'n_selected_fold_min': int(min(nsel_fold)), 'n_selected_fold_max': int(max(nsel_fold))}
    rB = {}
    yb = core.y_of(S[k]); assert np.array_equal(yb, y)
    for L in (1, 2, 3, 4):
        yf = fit_all(X[:, cols], y, L)
        qa, aa, _, _ = ao.alt_apparent(D1, G, L)
        qn, an, _ = ao.alt_nested_cv(D1, G, y, L)
        rA[str(L)] = dict(r2=q2(y, yf), auc_fit=roc_auc_score(y, yf), q2_app=qa, auc_app=aa, q2_nest=qn, auc_nest=an)
        yfb = fit_all(X12, y, L); yhb = loo_pls(X12, y, L)
        rB[str(L)] = dict(r2=q2(y, yfb), auc_fit=roc_auc_score(y, yfb), q2=q2(y, yhb), auc=roc_auc_score(y, yhb))
    out['A1'][k] = rA; out['B'][k] = rB
    print(k, 'A1 sel', len(cols), min(nsel_fold), max(nsel_fold),
          'R2', [round(rA[L]['r2'], 2) for L in '1234'], 'app', [round(rA[L]['q2_app'], 2) for L in '1234'],
          'nest', [round(rA[L]['q2_nest'], 2) for L in '1234'], '| B R2', [round(rB[L]['r2'], 2) for L in '1234'],
          'Q2', [round(rB[L]['q2'], 2) for L in '1234'], flush=True)
json.dump(out, open((_OUT + '/lv_routes.json'), 'w'), indent=1, default=float)
