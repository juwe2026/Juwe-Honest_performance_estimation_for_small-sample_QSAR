# -*- coding: utf-8 -*-
"""Second check of the one deviation found by mordred_check.py (benchmark/datasets/DS06, original rule, n = 30).
Recomputes the descriptor matrix with Mordred as mordred_check.py does and computes the observed values of every
arm of run_benchmark.py (hit counts, apparent and nested routes, controls, random forest) twice:
  (b) regenerated matrix with the archived correlation groups of DS06_info.json,
  (c) regenerated matrix with correlation groups regrouped from it by gruppieren() of datasets.py.
Both are compared with the archived observed values in ergebnisse/beobachtet.json, which were computed from the
archived matrix. Permutations are not recomputed.
Usage: python3 mordred_ds06_models.py <archive_root>/benchmark <archive_root>/benchmark/ergebnisse/beobachtet.json <out.json>
The first argument is the benchmark folder (or a copy holding routen.py, core_fast.py, run_benchmark.py, datasets.py and
datasets/DS06_*). The script reads from it and writes only the output file.
"""
import ast, json, os, sys, time
sys.dont_write_bytecode = True
from importlib.metadata import version
import numpy as np, pandas as pd
from scipy.stats import spearmanr
W, OBS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.chdir(W); sys.path.insert(0, W)
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
import run_benchmark as rb, routen
t0 = time.time()
tree = ast.parse(open('datasets.py', encoding='utf-8').read())
ns = {'np': np, 'spearmanr': spearmanr}
exec("from rdkit.Chem import Descriptors\nfrom mordred import Calculator, descriptors as mdesc\n"
     "CALC = Calculator(mdesc, ignore_3D=True)\nRD = [d for d in Descriptors._descList]\n", ns)
funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('deskriptoren', 'reinigen', 'gruppieren')]
exec(compile(ast.Module(body=funcs, type_ignores=[]), 'datasets.py', 'exec'), ns)
tab = pd.read_csv('datasets/DS06_y.csv'); X0 = pd.read_csv('datasets/DS06_X.csv')
X1, _ = ns['reinigen'](ns['deskriptoren']([Chem.MolFromSmiles(s) for s in tab.smiles]))
assert list(X1.columns) == list(X0.columns)
cols = list(X0.columns)
info = json.load(open('datasets/DS06_info.json'))
g_arch = [[cols.index(c) for c in grp if c in cols] for grp in info['gruppen']]
g_new, _ = ns['gruppieren'](X1)
y = tab['y'].values.astype(int)
A0 = X0.values.astype(float); A1 = X1.values.astype(float)
t_regen = time.time() - t0
from rdkit.Chem import rdFingerprintGenerator
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = np.array([np.array(gen.GetFingerprint(Chem.MolFromSmiles(s))) for s in tab['smiles']])

def alle(A, g):
    rb._D['DS06'] = (A, y, g, routen.zentralitaets_reps(A, g), fps, info['info'])
    return {art: rb.rechne('DS06', art, y) for art in rb.REIHENFOLGE}

def sel(A, g):
    """Selected descriptor names of the apparent routes on all compounds (same functions as the arms)."""
    repsf = routen.zentralitaets_reps(A, g)
    return dict(A1=sorted(cols[i] for i in rb.auswahl_a1(A, y, g)), A1u=sorted(cols[i] for i in rb.auswahl_a1u(A, y, repsf)),
                A2=sorted(cols[i] for i in rb.auswahl_a2(A, y)))

obs = json.load(open(OBS))['DS06']
res = dict(versions=dict(mordredcommunity=version('mordredcommunity'), rdkit=Chem.rdBase.rdkitVersion),
           n=len(y), p=len(cols), seconds_regeneration=round(t_regen, 1),
           max_relative_difference_matrix=float((np.abs(A1 - A0) / np.maximum(np.abs(A0), 1e-12)).max()),
           groups=dict(archived=len(g_arch), regrouped=len(g_new),
                       differing=len({frozenset(x) for x in g_arch} ^ {frozenset(x) for x in g_new})))
sel0 = sel(A0, g_arch)
for key, A, g in (('b_regenerated_archived_groups', A1, g_arch), ('c_regenerated_regrouped', A1, g_new)):
    ta = time.time(); r = alle(A, g)
    cmp = {}
    for art, v in r.items():
        o = obs[art]
        if art == 'Treffer':
            cmp[art] = dict(archived=o, recomputed=v, identical=o == v)
        else:
            cmp[art] = dict(archived=o, recomputed=list(v), abs_diff_q2=abs(v[0] - o[0]), abs_diff_auc=abs(v[1] - o[1]))
    s = sel(A, g)
    res[key] = dict(seconds=round(time.time() - ta, 1), arms=cmp,
                    apparent_selection_identical={k: s[k] == sel0[k] for k in s})
json.dump(res, open(OUT, 'w'), indent=1)
for key in ('b_regenerated_archived_groups', 'c_regenerated_regrouped'):
    print(key, res[key]['seconds'], 's', res[key]['apparent_selection_identical'])
    for art, c in res[key]['arms'].items():
        print('  ', art, c.get('identical', ''), c.get('abs_diff_q2', ''), c.get('abs_diff_auc', ''))
print('total', round(time.time() - t0, 1), 's')
