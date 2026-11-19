# -*- coding: utf-8 -*-
"""Effect of the exact duplicate columns on the correlation groups of the benchmark (main analysis, DS01 to DS16).
For every series: number of extra columns whose value vector equals that of an earlier column, whether every set of
identical columns lies in one archived correlation group (DSxx_info.json, key 'gruppen'), and the number of groups G
and the ratio G/p when gruppieren() of benchmark/datasets.py is applied to the archived matrix with and without the
duplicates. Both regroupings start from the CSV file read with float_precision='round_trip', which reproduces the archived groups
exactly, so the comparison isolates the duplicates.
Usage: python3 duplicate_groups_check.py <archive_root> <out.json>
"""
import ast, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT, OUT = sys.argv[1], sys.argv[2]
BM = os.path.join(ROOT, 'benchmark')
tree = ast.parse(open(os.path.join(BM, 'datasets.py'), encoding='utf-8').read())
ns = {'np': np, 'spearmanr': spearmanr}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'gruppieren'],
                        type_ignores=[]), 'datasets.py', 'exec'), ns)
res = {}
for k in range(1, 17):
    ds = 'DS%02d' % k; t0 = time.time()
    X = pd.read_csv(os.path.join(BM, 'datasets_k', ds + '_X.csv'), float_precision='round_trip')
    info = json.load(open(os.path.join(BM, 'datasets_k', ds + '_info.json')))
    A = X.values.astype(float); cols = list(X.columns)
    seen, sets = {}, {}
    for j in range(A.shape[1]):
        key = A[:, j].tobytes()
        if key in seen: sets.setdefault(seen[key], [seen[key]]).append(j)
        else: seen[key] = j
    extra = sorted(j for s in sets.values() for j in s[1:])
    gid = {c: i for i, g in enumerate(info['gruppen']) for c in g}
    same_group = all(len({gid.get(cols[j]) for j in s}) == 1 for s in sets.values())
    g_full, _ = ns['gruppieren'](X)
    Xd = X.drop(columns=[cols[j] for j in extra])
    g_dedup, _ = ns['gruppieren'](Xd)
    p, pd_ = X.shape[1], Xd.shape[1]
    res[ds] = dict(p=p, extra_duplicates=len(extra), duplicate_sets=len(sets),
                   every_set_in_one_archived_group=same_group, G_archived=len(info['gruppen']),
                   G_regrouped_csv=len(g_full), G_regrouped_without_duplicates=len(g_dedup),
                   Gp_archived=len(info['gruppen']) / p, Gp_regrouped_csv=len(g_full) / p,
                   Gp_without_duplicates=len(g_dedup) / pd_, seconds=round(time.time() - t0, 1))
    print(ds, res[ds], flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
