# -*- coding: utf-8 -*-
"""Regroup every archived benchmark matrix from its CSV file, read once with the default float
parser of pandas and once with float_precision='round_trip', and compare with the stored groups (DSxx_info.json).
Usage: python3 groups_csv_roundtrip.py <archive_root> <out.json>"""
import ast, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT, OUT = sys.argv[1:3]; BM = os.path.join(ROOT, 'benchmark')
tree = ast.parse(open(os.path.join(BM, 'datasets.py'), encoding='utf-8').read())
ns = {'np': np, 'spearmanr': spearmanr}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'gruppieren'],
                        type_ignores=[]), 'datasets.py', 'exec'), ns)
res = {}
for folder in ('datasets', 'datasets_k', 'datasets_m'):
    for f in sorted(os.listdir(os.path.join(BM, folder))):
        if not f.endswith('_X.csv'): continue
        key = '%s/%s' % (folder, f[:-6]); base = os.path.join(BM, folder, f[:-6])
        st = {frozenset(g) for g in json.load(open(base + '_info.json'))['gruppen']}
        r = {'n_groups_stored': len(st)}
        for name, fp in (('default', None), ('round_trip', 'round_trip')):
            X = pd.read_csv(base + '_X.csv', float_precision=fp); c = list(X.columns)
            g, _ = ns['gruppieren'](X); s = {frozenset(c[j] for j in gg) for gg in g}
            r['n_groups_' + name] = len(g); r[name + '_equals_stored'] = s == st
        res[key] = r; print(key, r, flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
print('default equal:', sum(r['default_equals_stored'] for r in res.values()), 'round_trip equal:',
      sum(r['round_trip_equals_stored'] for r in res.values()), 'of', len(res))
