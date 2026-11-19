"""Revision 2, Fig. 6a: apparent (not nested) Q2 of Routes A1, A1u, A2 with the code of Revision 1
(code/nested_generic.py, altorder.py, unsup_rep.py of this archive), so that the values are comparable with the
monoterpenoid values of Table 1.
Apparent = selection once on all compounds, then leave-one-out with this fixed selection
(for A1 the function altorder.alt_apparent, for A1u and A2 the same construction with unsup_select or
eff_filter + group_reps from the nested procedure).
Checks: (1) nested Q2 of the octapeptides against Additional file 2 (SecondDataSet_models),
(2) apparent and nested Q2 of the four evaluable monoterpenoid endpoints against Table 1.
Post hoc addition (Revision 2), not part of the pre-specification of the second data set (Additional file 5).
Usage: QSAR_ROUTEA_XLSX=<RouteA_matrix.xlsx> QSAR_RECALC_DIR=<tmp> python -I scheinbare_q2_rev1code.py <archive_root>/code <desc_pool.csv> <rows.json> <out.json>"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, sys.argv[1])
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
import nested_generic as ng, altorder as ao, unsup_rep as ur

def q2(y, yh): return float(1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2))
def loo_fixed(X, y, cols, nlv=1):
    N = len(y); yh = np.zeros(N)
    if len(cols) == 0: return np.full(N, 0.5)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xtr = X[tr][:, cols]; Xte = X[i:i+1][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(min(nlv, len(cols), int(tr.sum()) - 1), scale=False).fit((Xtr - mu) / sd, y[tr])
        yh[i] = m.predict((Xte - mu) / sd).ravel()[0]
    return yh
def apparent_all(D, G, REPS):
    X, y, N = D['X'], D['Y'], D['N']; R = rankdata(X, axis=0)
    a1 = ao.alt_apparent(D, G, 1)[0]
    c_u = ur.unsup_select(REPS, R, y, N)[0]; a1u = q2(y, loo_fixed(X, y, c_u))
    eff, p = ng.eff_filter(R, y, N); keep = np.where((eff > 0.3) & (p <= 0.05))[0]
    c2 = keep[ng.group_reps(X[:, keep], R[:, keep], eff[keep], N)] if len(keep) else keep
    a2 = q2(y, loo_fixed(X, y, c2))
    return dict(A1=a1, A1u=a1u, A2=a2), dict(A1=None, A1u=len(c_u), A2=len(c2))
def nested_all(D, G, REPS):
    y = D['Y']
    return dict(A1=float(ao.alt_nested_cv(D, G, y, 1)[0]), A1u=float(ur.unsup_nested_cv(D, REPS, y, 1)[0]),
                A2=float(ng.nested_cv(D, y, 1)[0]))
def make_D(X, y):
    N = len(y); F = []
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False; F.append(dict(tr=tr, i=i, Xtr=X[tr], R=rankdata(X[tr], axis=0), Xte=X[i:i+1]))
    return dict(X=X, Y=y, N=N, FOLDS=F)

out = {}
# ---- octapeptides, rule A (as revision/lone/run3.py)
Xp = pd.read_csv(sys.argv[2]); Xp = Xp[[c for c in Xp.columns if c != 'No']].values.astype(float)
rows = json.load(open(sys.argv[3])); THR = {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}
Gp, _ = ao.build_groups(Xp); REPSp = np.asarray(ur.centrality_reps(Xp, Gp))
ref_pep = {'EC': (0.1883, 0.0037, 0.1245), 'SA': (0.4461, 0.2515, 0.3941), 'PA': (0.2987, 0.3457, 0.2613), 'AB': (-0.2031, -0.3154, -0.2184)}
out['Octapeptide'] = {'gruppen': len(Gp)}
for k in ['EC', 'SA', 'PA', 'AB']:
    y = np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows])
    D = make_D(Xp, y); app, nsel = apparent_all(D, Gp, REPSp); nes = nested_all(D, Gp, REPSp)
    out['Octapeptide'][k] = {rt: dict(scheinbar=round(app[rt], 4), genestet=round(nes[rt], 4), genestet_AF2=ref_pep[k][i],
                                      optimismus=round(app[rt] - nes[rt], 4)) for i, rt in enumerate(['A1', 'A1u', 'A2'])}
    print('Pep', k, int(y.sum()), json.dumps(out['Octapeptide'][k]), flush=True)
# ---- monoterpenoids, four evaluable endpoints (check against Table 1)
YCOL = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPneu': 'S. Pneu (S. pneumoniae)', 'SPyo': 'S. pyo (S. pyogenes)'}
ref_app = {'HI': (0.40, 0.29, 0.46), 'SA': (0.51, 0.45, 0.52), 'SPneu': (0.56, 0.47, 0.55), 'SPyo': (0.49, 0.37, 0.53)}
ref_nes = {'HI': (0.0197, -0.2155, -0.1370), 'SA': (0.1042, 0.0811, 0.1550), 'SPneu': (0.3017, 0.1394, 0.1910), 'SPyo': (-0.1329, -0.2575, -0.1131)}
out['Monoterpenoide'] = {}
Gm = REPSm = None
for k, col in YCOL.items():
    D = ng.build('Cleaned_Steps1-4', ycol=col)
    if Gm is None: Gm = ao.get_groups(D['X']); REPSm = np.asarray(ur.centrality_reps(D['X'], Gm)); out['Monoterpenoide']['gruppen'] = len(Gm)
    app, nsel = apparent_all(D, Gm, REPSm); nes = nested_all(D, Gm, REPSm)
    out['Monoterpenoide'][k] = {rt: dict(scheinbar=round(app[rt], 4), scheinbar_Table1=ref_app[k][i], genestet=round(nes[rt], 4),
                                         genestet_Table1=ref_nes[k][i]) for i, rt in enumerate(['A1', 'A1u', 'A2'])}
    print('Mono', k, json.dumps(out['Monoterpenoide'][k]), flush=True)
json.dump(out, open(sys.argv[4], 'w'), indent=1)
