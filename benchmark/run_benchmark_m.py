"""SENSITIVITY ANALYSIS MEDIAN RULE of the main analysis (addendum, draft 25.10.2026): as run_benchmark_k.py, but datasets_m/, ergebnisse_m/, analysis codes + 300.

Originally: main run of the ChEMBL benchmark following the frozen specification.

Computes in packages and writes a file after every package. An interruption costs
at most the package in progress. Usage:  python3 run_benchmark_m.py <minutes>

Order: first all observed values, then the cheap permutations
(hit count, apparent models, controls, random forest), finally the
nested routes.
"""
import json, os, sys, time, glob
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
import routen
from core_fast import effekt_und_p
RDLogger.DisableLog('rdApp.*')

ERG = 'ergebnisse_m'
os.makedirs(ERG, exist_ok=True)
os.makedirs(f'{ERG}/perm', exist_ok=True)
CHUNK = 100
WORKER = 2
# sensitivity analysis median rule (main analysis): analysis codes + 100 (pre-specification) + 200 (supplementary specification) = + 300
CODE = {k: v + 300 for k, v in {'A1': 1, 'A1u': 2, 'A2': 3, 'ohneAuswahl': 4, 'ohneFilter': 5, 'RF': 6, 'Treffer': 7,
        'app_A1': 11, 'app_A1u': 12, 'app_A2': 13}.items()}
NPERM = {'A1': 2000, 'A1u': 2000, 'A2': 2000, 'ohneAuswahl': 2000, 'ohneFilter': 2000,
         'RF': 500, 'Treffer': 1000, 'app_A1': 1000, 'app_A1u': 1000, 'app_A2': 1000}
REIHENFOLGE = ['Treffer', 'app_A1', 'app_A1u', 'app_A2', 'ohneAuswahl', 'ohneFilter', 'RF', 'A1', 'A1u', 'A2']

_D = {}

def lade(ds):
    if ds in _D: return _D[ds]
    X = pd.read_csv(f'datasets_m/{ds}_X.csv'); cols = list(X.columns)
    A = X.values.astype(float)
    tab = pd.read_csv(f'datasets_m/{ds}_y.csv')
    y = tab['y'].values.astype(int)
    info = json.load(open(f'datasets_m/{ds}_info.json'))
    g = [[cols.index(c) for c in grp if c in cols] for grp in info['gruppen']]
    repsf = routen.zentralitaets_reps(A, g)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = np.array([np.array(gen.GetFingerprint(Chem.MolFromSmiles(s))) for s in tab['smiles']])
    _D[ds] = (A, y, g, repsf, fps, info['info'])
    return _D[ds]

def guete(y, pred):
    ss = ((y - y.mean()) ** 2).sum()
    q2 = 1 - ((y - pred) ** 2).sum() / ss if ss > 0 else 0.0
    try: auc = roc_auc_score(y, pred)
    except ValueError: auc = 0.5
    return float(q2), float(auc)

def pls_loocv(A, y, sel):
    n = len(y); pred = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        pred[i] = routen._pls_vorhersage(A[m], y[m], A[i], sel)
    return pred

def auswahl_a1(A, y, g):
    r, p = effekt_und_p(A, y)
    reps = np.array([grp[int(np.argmax(np.abs(r[grp])))] for grp in g])
    return reps[(p[reps] <= routen.P_MAX) & (np.abs(r[reps]) > routen.R_MIN)]

def auswahl_a1u(A, y, repsf):
    r, p = effekt_und_p(A, y)
    reps = np.asarray(repsf)
    return reps[(p[reps] <= routen.P_MAX) & (np.abs(r[reps]) > routen.R_MIN)]

def auswahl_a2(A, y):
    r, p = effekt_und_p(A, y)
    durch = np.where((p <= routen.P_MAX) & (np.abs(r) > routen.R_MIN))[0]
    if durch.size == 0: return durch
    gg = routen.gruppierung(A[:, durch])
    return np.array([durch[grp[int(np.argmax(np.abs(r[durch[grp]])))]] for grp in gg])

def treffer(A, y, g, repsf):
    r, p = effekt_und_p(A, y)
    reps1 = np.array([grp[int(np.argmax(np.abs(r[grp])))] for grp in g])
    reps2 = np.asarray(repsf)
    d = lambda idx: int(((p[idx] <= routen.P_MAX) & (np.abs(r[idx]) > routen.R_MIN)).sum())
    alle = int(((p <= routen.P_MAX) & (np.abs(r) > routen.R_MIN)).sum())
    return dict(A1=d(reps1), A1u=d(reps2), A2=alle)

def rechne(ds, art, y):
    A, _, g, repsf, fps, _ = lade(ds)
    if art == 'Treffer':  return treffer(A, y, g, repsf)
    if art == 'app_A1':   return guete(y, pls_loocv(A, y, auswahl_a1(A, y, g)))
    if art == 'app_A1u':  return guete(y, pls_loocv(A, y, auswahl_a1u(A, y, repsf)))
    if art == 'app_A2':   return guete(y, pls_loocv(A, y, auswahl_a2(A, y)))
    if art == 'ohneAuswahl': return guete(y, pls_loocv(A, y, np.arange(A.shape[1])))
    if art == 'ohneFilter':  return guete(y, pls_loocv(A, y, np.asarray(repsf)))
    if art == 'RF':
        rf = RandomForestClassifier(n_estimators=300, min_samples_leaf=1, oob_score=True,
                                    bootstrap=True, random_state=0, n_jobs=1)
        rf.fit(fps, y)
        pr = rf.oob_decision_function_[:, 1]
        pr = np.nan_to_num(pr, nan=0.5)
        return guete(y, pr)
    q2, pred = routen.loocv(A, y, g, 'A1' if art == 'A1' else ('A1u' if art == 'A1u' else 'A2'), repsf)
    return guete(y, pred)

def _job(arg):
    ds, art, dsnum, k = arg
    A, y0, *_ = lade(ds)
    rng = np.random.default_rng([93, dsnum, CODE[art], k])
    return rechne(ds, art, rng.permutation(y0))

def beobachtet(dss, nums):
    pfad = f'{ERG}/beobachtet.json'
    out = json.load(open(pfad)) if os.path.exists(pfad) else {}
    neu = False
    for ds in dss:
        out.setdefault(ds, {})
        A, y, *_ = lade(ds)
        for art in REIHENFOLGE:
            if art in out[ds]: continue
            out[ds][art] = rechne(ds, art, y); neu = True
    if neu: json.dump(out, open(pfad, 'w'), indent=1)
    return out

if __name__ == '__main__':
    minuten = float(sys.argv[1]) if len(sys.argv) > 1 else 110
    ende = time.time() + minuten * 60
    dss = sorted({f.split('/')[-1][:4] for f in glob.glob('datasets_m/DS*_X.csv')})
    nums = {ds: int(ds[2:]) for ds in dss}
    t0 = time.time()
    beobachtet(dss, nums)
    print(f'observed values complete ({time.time()-t0:.0f}s)', flush=True)
    pakete = [(ds, art, s) for art in REIHENFOLGE for ds in dss
              for s in range(0, NPERM[art], CHUNK)]
    offen = [p for p in pakete if not os.path.exists(f'{ERG}/perm/{p[0]}_{p[1]}_{p[2]:04d}.json')]
    print(f'packages total {len(pakete)}, open {len(offen)}', flush=True)
    fertig = 0
    with ProcessPoolExecutor(max_workers=WORKER) as ex:
        for ds, art, s in offen:
            if time.time() > ende:
                print('time limit reached', flush=True); break
            ta = time.time()
            args = [(ds, art, nums[ds], k) for k in range(s, min(s + CHUNK, NPERM[art]))]
            werte = list(ex.map(_job, args))
            json.dump(dict(ds=ds, art=art, start=s, n=len(werte), werte=werte),
                      open(f'{ERG}/perm/{ds}_{art}_{s:04d}.json', 'w'))
            fertig += 1
            print(f'{ds} {art} {s:04d}  {len(werte)} permutations  {time.time()-ta:.0f}s', flush=True)
    print(f'DONE in this block: {fertig} packages, open afterwards: {len(offen)-fertig}', flush=True)
