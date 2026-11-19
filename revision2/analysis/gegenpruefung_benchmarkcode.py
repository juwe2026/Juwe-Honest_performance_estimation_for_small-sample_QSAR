"""Revision 2, cross-check: the three data sets of the study (monoterpenoids, four evaluable endpoints,
octapeptides rule A, Gaucher) with the CODE OF THE BENCHMARK (routen.py, core_fast.py, selection functions from
run_benchmark_k.py) instead of the code of Revision 1. Purpose: check how much the differences between the implementations
(tie correction in the Mann-Whitney p, choice of representative under A2, centrality measure under A1u) shift the
calibrated inflation factors and the optimism, so that Fig. 6 compares like with like.
Null distribution of the hit count: 1,000 permutations, seeds default_rng([93, 900 + data set, 7, k]).
Calibrated factor as in the benchmark (auswertung.py): SD(null, ddof=1) / sqrt(m q (1-q)), q = mean(null)/m,
m = number of descriptors (A2) or groups (A1, A1u).
Usage: QSAR_ROUTEA_XLSX=<RouteA_matrix.xlsx> QSAR_RECALC_DIR=<tmp> python -I gegenpruefung_benchmarkcode.py <benchmark_scripts> <archive_root> <out.json>"""
import sys, os, json, time, numpy as np, pandas as pd
BM, AF7, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, BM)
import routen
from core_fast import effekt_und_p
sys.path.insert(1, os.path.join(AF7, 'code'))
import nested_generic as ng, gaucher_pipe as gp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gruppe_test import gruppieren_vec   # = datasets.gruppieren of the benchmark (checked: identical)

def auswahl_a1(A, y, g):
    r, p = effekt_und_p(A, y); reps = np.array([grp[int(np.argmax(np.abs(r[grp])))] for grp in g])
    return reps[(p[reps] <= routen.P_MAX) & (np.abs(r[reps]) > routen.R_MIN)]
def auswahl_a1u(A, y, repsf):
    r, p = effekt_und_p(A, y); reps = np.asarray(repsf)
    return reps[(p[reps] <= routen.P_MAX) & (np.abs(r[reps]) > routen.R_MIN)]
def auswahl_a2(A, y):
    r, p = effekt_und_p(A, y); durch = np.where((p <= routen.P_MAX) & (np.abs(r) > routen.R_MIN))[0]
    if durch.size == 0: return durch
    gg = routen.gruppierung(A[:, durch])
    return np.array([durch[grp[int(np.argmax(np.abs(r[durch[grp]])))]] for grp in gg])
def pls_loocv(A, y, sel):
    n = len(y); pred = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False; pred[i] = routen._pls_vorhersage(A[m], y[m], A[i], sel)
    return pred
def q2(y, pr): return float(1 - ((y - pr) ** 2).sum() / ((y - y.mean()) ** 2).sum())
def treffer(A, y, g, repsf):
    r, p = effekt_und_p(A, y)
    reps1 = np.array([grp[int(np.argmax(np.abs(r[grp])))] for grp in g]); reps2 = np.asarray(repsf)
    d = lambda idx: int(((p[idx] <= routen.P_MAX) & (np.abs(r[idx]) > routen.R_MIN)).sum())
    return dict(A1=d(reps1), A1u=d(reps2), A2=int(((p <= routen.P_MAX) & (np.abs(r) > routen.R_MIN)).sum()))

def lauf(name, dsid, X, endpunkte, nperm=1000):
    t0 = time.time()
    g = gruppieren_vec(X)                           # label-free groups as in the benchmark (datasets.gruppieren)
    repsf = routen.zentralitaets_reps(X, g)
    G, p = len(g), X.shape[1]; res = dict(G=G, p=p, endpunkte={})
    for ep, y in endpunkte.items():
        obs = treffer(X, y, g, repsf)
        null = {r: [] for r in ['A1', 'A1u', 'A2']}
        for k in range(nperm):
            t = treffer(X, np.random.default_rng([93, 900 + dsid, 7, k]).permutation(y), g, repsf)
            for r in null: null[r].append(t[r])
        fak = {}
        for r, m in [('A1', G), ('A1u', G), ('A2', p)]:
            nl = np.array(null[r], float); qn = nl.mean() / m; sd = nl.std(ddof=1)
            fak[r] = dict(kalibriert=round(float(sd / np.sqrt(m * qn * (1 - qn))), 3), nominal=round(float(sd / np.sqrt(m * 0.05 * 0.95)), 3),
                          null_rate=round(float(qn), 4), beobachtet=obs[r])
        q = {}
        for r, sel in [('A1', auswahl_a1(X, y, g)), ('A1u', auswahl_a1u(X, y, repsf)), ('A2', auswahl_a2(X, y))]:
            app = q2(y, pls_loocv(X, y, sel)); nes = q2(y, routen.loocv(X, y, g, r, repsf)[1])
            q[r] = dict(scheinbar=round(app, 4), genestet=round(nes, 4), optimismus=round(app - nes, 4))
        res['endpunkte'][ep] = dict(aktiv=int(y.sum()), faktoren=fak, q2=q)
        print(name, ep, 'G', G, json.dumps(fak), json.dumps(q), f'{time.time()-t0:.0f}s', flush=True)
    return res

out = {}
# octapeptides, rule A
Xp = pd.read_csv(os.path.join(AF7, 'revision/lone/desc_pool.csv')); Xp = Xp[[c for c in Xp.columns if c != 'No']].values.astype(float)
rows = json.load(open(os.path.join(AF7, 'revision/lone/rows.json'))); THR = {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}
ep = {k: np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows]) for k in THR}
out['Octapeptide'] = lauf('Octapeptide', 2, Xp, ep)
# Gaucher
X, y, mz, _ = gp.load(os.path.join(AF7, 'data/Gaucherdata_paper39_median_normalised.csv')); Xg, _, _ = gp.clean(X, mz)
out['Gaucher'] = lauf('Gaucher', 3, Xg, {'Gaucher': y})
# monoterpenoids, four evaluable endpoints
YCOL = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPneu': 'S. Pneu (S. pneumoniae)', 'SPyo': 'S. pyo (S. pyogenes)'}
Ds = {k: ng.build('Cleaned_Steps1-4', ycol=c) for k, c in YCOL.items()}
Xm = Ds['HI']['X']
assert all(np.array_equal(Ds[k]['X'], Xm, equal_nan=True) for k in Ds)
out['Monoterpenoide'] = lauf('Monoterpenoide', 1, Xm, {k: Ds[k]['Y'] for k in Ds})
json.dump(out, open(OUT, 'w'), indent=1)
print('done')
