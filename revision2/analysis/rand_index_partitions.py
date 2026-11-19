# -*- coding: utf-8 -*-
"""Rand index and adjusted Rand index of the fold-internal correlation partitions of the monoterpenoid pool against
the global partition (Additional file 1, section before Table S2, item G2 of 16.11.2026).

The global partition groups the 2,208 cleaned descriptors (Additional file 3) into 182 correlation groups with
code/altorder.py. For each of the 31 leave-one-out folds the partition is rebuilt from the 30 training compounds and
compared with the global one. The script also gives the Rand index of the trivial all-singletons partition, which
shows how much of the Rand index comes from the large majority of descriptor pairs that fall between groups.
The adjusted Rand index (Hubert and Arabie) corrects for this chance agreement. Both are taken from scikit-learn.

Usage:  python rand_index_partitions.py <archive_root> <empty_work_dir> <out.json>
Additional file 3 must be placed in data/ as described in the README. Nothing is written into the archive.
Verified with Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, scikit-learn 1.9.1. Running time about two minutes.
"""
import sys, os, json, time
R = sys.argv[1]; os.environ['QSAR_RECALC_DIR'] = sys.argv[2]
os.environ['QSAR_ROUTEA_XLSX'] = R + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'
sys.path.insert(0, R + '/code')
import numpy as np
from sklearn.metrics import adjusted_rand_score, rand_score
import nested_generic as ng, altorder as ao

D = ng.build('Cleaned_Steps1-4'); X = D['X']; N, P = X.shape


def labels(groups):
    g = np.empty(P, dtype=np.int64)
    for k, members in enumerate(groups): g[np.asarray(members, dtype=int)] = k
    return g


t0 = time.time()
G, _ = ao.build_groups(X); gl = labels(G)
pairs = P * (P - 1) // 2; within = int(sum(len(g) * (len(g) - 1) // 2 for g in G))
single = (float(rand_score(gl, np.arange(P))), float(adjusted_rand_score(gl, np.arange(P))))
print('global partition: %d groups of %d descriptors, %d of %d pairs within groups' % (len(G), P, within, pairs))
print('all-singletons partition against the global one: Rand index %.6f, adjusted %.6f' % single)
folds = []
for i in range(N):
    tr = np.ones(N, bool); tr[i] = False
    g, _ = ao.build_groups(X[tr]); gf = labels(g)
    folds.append(dict(left_out=i, groups=len(g), rand=float(rand_score(gl, gf)), adjusted_rand=float(adjusted_rand_score(gl, gf))))
    print('fold %2d: %3d groups, Rand index %.4f, adjusted %.4f (%.0f s)' % (i, len(g), folds[-1]['rand'], folds[-1]['adjusted_rand'], time.time() - t0), flush=True)
r = np.array([f['rand'] for f in folds]); a = np.array([f['adjusted_rand'] for f in folds]); ng_ = [f['groups'] for f in folds]
summary = dict(groups_global=len(G), descriptors=P, pairs=pairs, pairs_within_groups=within,
               rand_all_singletons=single[0], adjusted_rand_all_singletons=single[1],
               groups_folds_min=min(ng_), groups_folds_max=max(ng_), rand_min=float(r.min()), rand_max=float(r.max()),
               adjusted_rand_min=float(a.min()), adjusted_rand_max=float(a.max()), adjusted_rand_median=float(np.median(a)),
               adjusted_rand_mean=float(a.mean()), folds_adjusted_rand_below_0_6=int((a < 0.6).sum()))
print(json.dumps(summary, indent=1))
json.dump(dict(summary=summary, folds=folds), open(sys.argv[3], 'w'), indent=1)
