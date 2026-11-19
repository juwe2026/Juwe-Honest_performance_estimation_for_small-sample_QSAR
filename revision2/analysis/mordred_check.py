# -*- coding: utf-8 -*-
"""Check of the Mordred version for the benchmark descriptor matrices (6 November 2026).

The version of Mordred that computed the descriptor matrices of the benchmark (benchmark/datasets.py) was not
recorded. This script recomputes every matrix from the archived structures (column smiles of DSxx_y.csv) with the
functions deskriptoren() and reinigen() of benchmark/datasets.py, taken from its source unchanged, and compares the
result with the archived DSxx_X.csv: number of raw descriptors, cleaned columns (names and order) and values.
The structures are re-read from the canonical SMILES, so the atom order can differ from the ChEMBL records.
Usage: python3 mordred_check.py <archive_root> <out.json>
"""
import ast, json, os, sys, time
from importlib.metadata import version
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
ROOT, OUT = sys.argv[1], sys.argv[2]
BM = os.path.join(ROOT, 'benchmark')
tree = ast.parse(open(os.path.join(BM, 'datasets.py'), encoding='utf-8').read())
ns = {}
exec("from rdkit.Chem import Descriptors\nfrom mordred import Calculator, descriptors as mdesc\n"
     "CALC = Calculator(mdesc, ignore_3D=True)\nRD = [d for d in Descriptors._descList]\n", ns)
funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('deskriptoren', 'reinigen')]
exec(compile(ast.Module(body=funcs, type_ignores=[]), 'datasets.py', 'exec'), ns)
out = dict(versions=dict(mordredcommunity=version('mordredcommunity'), rdkit=Chem.rdBase.rdkitVersion,
                         numpy=np.__version__, pandas=pd.__version__), data_sets={})
for folder in ('datasets', 'datasets_k', 'datasets_m'):
    for f in sorted(os.listdir(os.path.join(BM, folder))):
        if not f.endswith('_X.csv'): continue
        ds = f[:-6]; t0 = time.time()
        y = pd.read_csv(os.path.join(BM, folder, ds + '_y.csv')); X0 = pd.read_csv(os.path.join(BM, folder, f))
        info = json.load(open(os.path.join(BM, folder, ds + '_info.json')))['info']
        X, raw = ns['reinigen'](ns['deskriptoren']([Chem.MolFromSmiles(s) for s in y.smiles]))
        r = dict(n=len(y), raw=[raw, info['p_roh']], cleaned=[X.shape[1], X0.shape[1]], same_columns=list(X.columns) == list(X0.columns))
        if r['same_columns']:
            A, B = X.values.astype(float), X0.values.astype(float)
            rel = np.abs(A - B) / np.maximum(np.abs(B), 1e-12)
            r['share_identical_cells'] = round(float((A == B).mean()), 4)
            r['max_relative_difference'] = float(rel.max())
            r['descriptors_above_1e-9'] = sorted({X.columns[j] for i, j in np.argwhere(rel > 1e-9)})
        out['data_sets']['%s/%s' % (folder, ds)] = r
        print(folder, ds, r, '%.1fs' % (time.time() - t0), flush=True)
json.dump(out, open(OUT, 'w'), indent=1)
