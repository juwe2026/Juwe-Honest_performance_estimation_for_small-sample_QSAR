# -*- coding: utf-8 -*-
"""Descriptive sensitivity analysis. Hit-count null distributions, calibrated inflation factors and the tests of
H2 and H3 for the 16 natural series of the main analysis, once as archived and once without the exact duplicate
columns. Not pre-specified, descriptive only.

Step 1 (check of the reimplementation): with the archived matrices, groups and permutation seeds the observed hit
counts and all 1,000 null hit counts per series must equal ergebnisse_k/beobachtet.json and ergebnisse_k/perm/.
Step 2: every column whose value vector equals that of an earlier column is dropped. The correlation groups are
regrouped with gruppieren() of datasets.py from the CSV file read with float_precision='round_trip' (this reproduces
the archived groups exactly when nothing is dropped), the label-free representatives are recomputed with
routen.zentralitaets_reps(), and the same permutation streams are used. The factors follow auswertung_k.py.
Usage: python3 duplicate_hits_check.py <archive_root> <out.json>
"""
import ast, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr, page_trend_test, wilcoxon
ROOT, OUT = sys.argv[1:3]
BM = os.path.join(ROOT, 'benchmark'); os.chdir(BM); sys.path.insert(0, BM)
import routen
from core_fast import effekt_und_p
tree = ast.parse(open('datasets.py', encoding='utf-8').read())
ns = {'np': np, 'spearmanr': spearmanr}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'gruppieren'],
                        type_ignores=[]), 'datasets.py', 'exec'), ns)
src = open('auswertung_k.py', encoding='utf-8').read()
PRIMAER = ['DS%02d' % k for k in range(1, 17)]
CODE_T = 207  # 'Treffer' + 200, run_benchmark_k.py
NPERM = 1000
ROUTEN = ['A1', 'A1u', 'A2']


def treffer(A, y, g, repsf):
    """Identical with treffer() in run_benchmark.py."""
    r, p = effekt_und_p(A, y)
    reps1 = np.array([grp[int(np.argmax(np.abs(r[grp])))] for grp in g])
    reps2 = np.asarray(repsf)
    d = lambda idx: int(((p[idx] <= routen.P_MAX) & (np.abs(r[idx]) > routen.R_MIN)).sum())
    alle = int(((p <= routen.P_MAX) & (np.abs(r) > routen.R_MIN)).sum())
    return dict(A1=d(reps1), A1u=d(reps2), A2=alle)


def faktoren(obs, nulls, tests):
    out = {}
    for r in ROUTEN:
        nl = np.array([w[r] for w in nulls], float); m_ = tests[r]
        sd = np.std(nl, ddof=1); qn = nl.mean() / m_
        out[r] = dict(beobachtet=obs[r], null_mittel=float(nl.mean()), null_sd=float(sd), tests=m_,
                      inf_kalibriert=float(sd / np.sqrt(m_ * qn * (1 - qn))) if 0 < qn < 1 else float('nan'),
                      p=float((1 + np.sum(nl >= obs[r])) / (1 + len(nl))))
    return out


def archived_null(ds):
    v = []
    for s in range(0, NPERM, 100):
        d = json.load(open('ergebnisse_k/perm/%s_Treffer_%04d.json' % (ds, s)))
        assert d['start'] == s; v += d['werte']
    return v


beob = json.load(open('ergebnisse_k/beobachtet.json'))
ausw = json.load(open('ergebnisse_k/auswertung.json'))
res = dict(series={}, note='descriptive sensitivity analysis, not pre-specified')
for ds in PRIMAER:
    t0 = time.time(); dsnum = int(ds[2:])
    X = pd.read_csv('datasets_k/%s_X.csv' % ds); cols = list(X.columns); A = X.values.astype(float)
    y = pd.read_csv('datasets_k/%s_y.csv' % ds)['y'].values.astype(int)
    info = json.load(open('datasets_k/%s_info.json' % ds))
    g = [[cols.index(c) for c in grp if c in cols] for grp in info['gruppen']]
    repsf = routen.zentralitaets_reps(A, g)
    arch = archived_null(ds)
    # The main analysis took four series unchanged from the original rule (README_benchmark.md), whose permutation
    # streams use the analysis code without the offset of 200. The code is chosen by the first 100 archived draws.
    code = None
    for c in (CODE_T, CODE_T - 200):
        pr = [np.random.default_rng([93, dsnum, c, k]).permutation(y) for k in range(100)]
        if [treffer(A, yv, g, repsf) for yv in pr] == arch[:100]: code = c; break
    assert code is not None, ds
    perms = [np.random.default_rng([93, dsnum, code, k]).permutation(y) for k in range(NPERM)]
    # step 1
    obs0 = treffer(A, y, g, repsf); null0 = [treffer(A, yv, g, repsf) for yv in perms]
    check = dict(seed_code=code, observed_equal=obs0 == beob[ds]['Treffer'], null_equal=null0 == arch)
    tests0 = dict(A1=info['info']['gruppen'], A1u=info['info']['gruppen'], A2=info['info']['p_sauber'])
    f0 = faktoren(obs0, null0, tests0)
    check['factors_equal'] = all(abs(f0[r]['inf_kalibriert'] - ausw['tabelle'][ds]['T_' + r]['inf_kalibriert']) < 1e-12
                                 for r in ROUTEN)
    # step 2
    seen, drop = {}, []
    for j in range(A.shape[1]):
        key = A[:, j].tobytes()
        if key in seen: drop.append(j)
        else: seen[key] = j
    keep = [j for j in range(A.shape[1]) if j not in set(drop)]
    Ad = A[:, keep]
    Xr = pd.read_csv('datasets_k/%s_X.csv' % ds, float_precision='round_trip').iloc[:, keep]
    gd, _ = ns['gruppieren'](Xr)
    gd = [list(map(int, grp)) for grp in gd]
    repsfd = routen.zentralitaets_reps(Ad, gd)
    obs1 = treffer(Ad, y, gd, repsfd); null1 = [treffer(Ad, yv, gd, repsfd) for yv in perms]
    tests1 = dict(A1=len(gd), A1u=len(gd), A2=Ad.shape[1])
    f1 = faktoren(obs1, null1, tests1)
    res['series'][ds] = dict(check=check, p=A.shape[1], duplicates=len(drop), G=len(g), p_without=Ad.shape[1],
                             G_without=len(gd), Gp=info['info']['redundanz_G_p'], Gp_without=len(gd) / Ad.shape[1],
                             archived=f0, without_duplicates=f1, seconds=round(time.time() - t0, 1))
    print(ds, check, 'IF A2 %.3f -> %.3f' % (f0['A2']['inf_kalibriert'], f1['A2']['inf_kalibriert']),
          'A1 %.3f -> %.3f' % (f0['A1']['inf_kalibriert'], f1['A1']['inf_kalibriert']),
          'A1u %.3f -> %.3f' % (f0['A1u']['inf_kalibriert'], f1['A1u']['inf_kalibriert']),
          '%.1fs' % (time.time() - t0), flush=True)

S = res['series']
for tag, key, gpkey in (('archived', 'archived', 'Gp'), ('without_duplicates', 'without_duplicates', 'Gp_without')):
    IF = {r: np.array([S[d][key][r]['inf_kalibriert'] for d in PRIMAER]) for r in ROUTEN}
    pt = page_trend_test(np.column_stack([IF['A1u'], IF['A1'], IF['A2']]))
    Gp = np.array([S[d][gpkey] for d in PRIMAER])
    sp = {r: spearmanr(Gp, np.log(IF[r]), alternative='less') for r in ROUTEN}
    res[tag] = dict(H2=dict(L=float(pt.statistic), p=float(pt.pvalue),
                            A2_gt_A1=float(wilcoxon(IF['A2'] - IF['A1'], alternative='greater').pvalue),
                            A1_gt_A1u=float(wilcoxon(IF['A1'] - IF['A1u'], alternative='greater').pvalue),
                            order_complete=int(((IF['A2'] > IF['A1']) & (IF['A1'] > IF['A1u'])).sum()),
                            median_IF={r: float(np.median(IF[r])) for r in ROUTEN}),
                    H3=dict(rho=float(sp['A2'].statistic), p=float(sp['A2'].pvalue),
                            secondary={r: dict(rho=float(sp[r].statistic), p=float(sp[r].pvalue)) for r in ['A1', 'A1u']}))
H = ausw['hypothesen']
res['archived_check'] = dict(H2_p_equal=abs(res['archived']['H2']['p'] - H['H2']['p']) < 1e-12,
                             H3_rho_equal=abs(res['archived']['H3']['rho'] - H['H3']['rho']) < 1e-12)
json.dump(res, open(OUT, 'w'), indent=1)
for tag in ('archived', 'without_duplicates'):
    print(tag, json.dumps(res[tag]))
print('archived check', res['archived_check'])
