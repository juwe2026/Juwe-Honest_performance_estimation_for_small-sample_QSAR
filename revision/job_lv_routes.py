# Revision (EW 03.10.2026): Table S5 (continued) fuer Route A1 und Route B, 1-4 LV.
# R2 = Anpassung an alle 31 Verbindungen (Route A1: Auswahl auf allen Verbindungen), wie Table S10;
# Q2 apparent = Auswahl einmal auf allen Verbindungen, danach LOO (Route A1, wie Table 1);
# Q2 nested = Auswahl in jedem Fold wiederholt (Route A1); Route B: feste 12 Deskriptoren, LOO.
# Deterministisch, keine Permutationen. Gleicher Code wie jobI2_fallback.py (fixcode).
import sys, json
sys.path.insert(0, '/home/claude/rev/fixcode')
exec(open('/home/claude/rev/common.py').read())
for _m in ['nested_generic', 'altorder', 'unsup_rep', 'core']: sys.modules.pop(_m, None)
sys.path.insert(0, '/home/claude/rev/fixcode')
import altorder as ao
assert '/home/claude/rev/fixcode' in ao.__file__
sys.path.insert(1, '/home/claude/rev/QSAR_selection_leakage/code')
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
json.dump(out, open('/home/claude/rev/out/lv_routes.json', 'w'), indent=1, default=float)
