# -*- coding: utf-8 -*-
"""Effect of the one deviation found by mordred_check.py: data set DS06 of the analysis under the original rule
(benchmark/datasets/DS06). Mordred 2.0.7 reproduces this matrix with a maximum relative difference of 0.2 % in five
descriptors. This script recomputes the matrix as mordred_check.py does and checks whether the difference reaches
the label-free steps of the pipeline: the ranks of every descriptor (all filters and the univariate tests are
rank-based), the Pearson and Spearman correlation graph at |r| > 0.7 and the correlation groups of gruppieren() in
benchmark/datasets.py. Usage: python3 mordred_ds06_check.py <archive_root> <out.json>
"""
import ast, json, os, sys
from importlib.metadata import version
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
ROOT, OUT = sys.argv[1], sys.argv[2]
BM = os.path.join(ROOT, 'benchmark')
tree = ast.parse(open(os.path.join(BM, 'datasets.py'), encoding='utf-8').read())
ns = {'np': np, 'spearmanr': spearmanr}
exec("from rdkit.Chem import Descriptors\nfrom mordred import Calculator, descriptors as mdesc\n"
     "CALC = Calculator(mdesc, ignore_3D=True)\nRD = [d for d in Descriptors._descList]\n", ns)
funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('deskriptoren', 'reinigen', 'gruppieren')]
exec(compile(ast.Module(body=funcs, type_ignores=[]), 'datasets.py', 'exec'), ns)
y = pd.read_csv(os.path.join(BM, 'datasets', 'DS06_y.csv')); X0 = pd.read_csv(os.path.join(BM, 'datasets', 'DS06_X.csv'))
X, raw = ns['reinigen'](ns['deskriptoren']([Chem.MolFromSmiles(s) for s in y.smiles]))
assert list(X.columns) == list(X0.columns)
A, B = X.values.astype(float), X0.values.astype(float)
rel = np.abs(A - B) / np.maximum(np.abs(B), 1e-12)
cols = sorted({X.columns[j] for i, j in np.argwhere(rel > 1e-9)})
res = dict(versions=dict(mordredcommunity=version('mordredcommunity'), rdkit=Chem.rdBase.rdkitVersion),
           n=len(y), p=X.shape[1], max_relative_difference=float(rel.max()), descriptors=cols,
           compounds_affected=int((rel > 1e-9).any(axis=1).sum()))
res['ranks_identical'] = {c: bool(np.array_equal(rankdata(X[c]), rankdata(X0[c]))) for c in cols}
P0, P1 = np.corrcoef(B, rowvar=False), np.corrcoef(A, rowvar=False)
S0, S1 = spearmanr(B).statistic, spearmanr(A).statistic
G0 = (np.abs(np.nan_to_num(P0)) > 0.7) | (np.abs(np.nan_to_num(S0)) > 0.7)
G1 = (np.abs(np.nan_to_num(P1)) > 0.7) | (np.abs(np.nan_to_num(S1)) > 0.7)
res['correlation_graph_identical'] = bool(np.array_equal(G0, G1))
res['max_abs_pearson_change'] = float(np.nanmax(np.abs(P1 - P0)))
g0, g1 = ns['gruppieren'](X0), ns['gruppieren'](X)
(gr0, c0), (gr1, c1) = g0, g1
res['groups_identical'] = gr0 == gr1
d = np.argwhere(np.triu(G0 != G1, 1))
res['edges_changed'] = [[X.columns[i], X.columns[j], round(float(P0[i, j]), 6), round(float(P1[i, j]), 6),
                         round(float(S0[i, j]), 6), round(float(S1[i, j]), 6)] for i, j in d]
res['n_groups'] = [len(gr0), len(gr1)]
s0, s1 = {frozenset(g) for g in gr0}, {frozenset(g) for g in gr1}
res['groups_differing'] = [len(s0 - s1), len(s1 - s0)]
res['hubs_differing'] = sorted({X.columns[g[0]] for g in gr0} ^ {X.columns[g[0]] for g in gr1})
info = json.load(open(os.path.join(BM, 'datasets', 'DS06_info.json')))['info']
res['archived_info_groups'] = {k: v for k, v in info.items() if 'grupp' in k.lower() or k in ('G', 'g', 'p', 'n')}
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
