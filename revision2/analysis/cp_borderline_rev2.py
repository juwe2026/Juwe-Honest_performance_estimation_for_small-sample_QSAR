# -*- coding: utf-8 -*-
"""Clopper-Pearson intervals (95 %) for the Monte Carlo uncertainty of the reported permutation p-values (Table S10).

Same items and the same method as cp_borderline.py. Changed are the thresholds, the paths and one omission:
  - thresholds 0.05 and 0.05/3 = 0.0167, the Bonferroni threshold for the three effective endpoints of each
    MIC series (Table S19). The former threshold 0.0125 (four evaluable endpoints) is no longer used.
    For transparency the output still flags the items whose interval contains 0.0125.
  - paths relative to the root of the code archive, so the script runs from any directory.
  - Table 5 now includes P. aeruginosa of the octapeptides (14 actives of 37, an evaluable endpoint). cp_borderline.py
    looped over E. coli, S. aureus and A. baumannii only, so this endpoint was missing from the earlier Table S10.

k = number of permuted values >= observed (ties counted, tolerance 1e-9), n = number of permutations.
The reported p is (1 + k)/(1 + n). The interval is the exact binomial interval for k/n, the proportion
that the Monte Carlo procedure estimates.

Usage:  python revision2/analysis/cp_borderline_rev2.py [archive_root] [out.json]
"""
import json, os, sys
import numpy as np
from scipy.stats import beta

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'results', 'revision', 'cp_borderline_rev2.json')
O = os.path.join(ROOT, 'results', 'revision') + os.sep
RP = os.path.join(ROOT, 'results', 'raw_permutations') + os.sep
LONE = os.path.join(O, 'lone') + os.sep
J = lambda f: json.load(open(f))
THR = (0.05, 0.05 / 3)          # decision thresholds used in the manuscript
THR_OLD = 0.0125                 # reported only for the comparison with the earlier table

def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)

def k_of(p, n):
    k = p * (n + 1) - 1
    assert abs(k - round(k)) < 1e-6, (p, n)
    return int(round(k))

def arch(prefix, s, m):
    z = np.load(RP + '%s_%s.npz' % (prefix, s))
    nul, ob = (z['q2s'], float(z['obs_q2'])) if m == 'q2' else (z['aucs'], float(z['obs_auc']))
    return float((1 + np.sum(nul >= ob - 1e-9)) / (1 + len(nul))), len(nul)

SN = {'HI': 'H. influenzae', 'SA': 'S. aureus', 'SPneu': 'S. pneumoniae', 'SPyo': 'S. pyogenes',
      'EC': 'E. coli', 'AB': 'A. baumannii', 'PA': 'P. aeruginosa'}   # PA only for the octapeptides (Table 5)
ITEMS = []
def add(src, model, s, m, p, n):
    ITEMS.append(dict(Source=src, Model=model, Endpoint=SN[s], Metric='Q²' if m == 'q2' else 'AUC', p=p, n=n))

T1 = J(O + 'I_fallback_seeds.json'); A = J(O + 'A_nofilter.json'); RF = J(O + 'D_rf.json')
for s in ('HI', 'SA', 'SPneu', 'SPyo'):
    for m in ('q2', 'auc'):
        for r in ('A1', 'A1u', 'A2'):
            add('Table 1', 'Route %s' % r, s, m, T1[s][r]['p_' + m], 2000)
        add('Table 1', 'No selection, all 2,208 descriptors', s, m, A['all2208'][s]['p_' + m], 2000)
        add('Table 1', 'No filter, 182 label-free representatives', s, m, A['reps182_centrality'][s]['p_' + m], 2000)
        add('Table 1', 'Random forest, out-of-bag', s, m, RF['RouteA_all2208'][s]['p_' + m], 500)
        for pre, lab in (('permB_12', 'Route B, 12 descriptors, 1 LV'), ('permB_12nlv3', 'Route B, 12 descriptors, 3 LV')):
            p, n = arch(pre, s, m); add('Table 1', lab, s, m, p, n)
        _b = J(O + 'B25_cineole0_seed701.json')[s]; add('Table 2', '25 descriptors', s, m, _b['p_' + m], 2000)
        add('Table 3', 'Route A1 top-15 by effect size, nested', s, m, J(O + 'cache_topk_seeds/%s_A1_top15.json' % s)['p_' + m], 2000)
        add('Table 4', 'Route A1 top-15 by VIP, nested', s, m, J(O + 'cache_topk_seeds/%s_A1_viptop1.json' % s)['p_' + m], 2000)
for rule, lab in (('A', ''), ('C', ', sensitivity definition')):
    R = J(LONE + 'res_%s.json' % rule); RX = J(LONE + 'res_%s_fix.json' % rule)
    for s in ('EC', 'SA', 'PA', 'AB'):
        for m in ('q2', 'auc'):
            for r in ('A1', 'A1u', 'A2'):
                add('Table 5', 'Route %s%s' % (r, lab), s, m, RX[s][r]['p_' + m], 2000)
            add('Table 5', 'No selection, all 1,275 descriptors%s' % lab, s, m, R[s]['noSelection_all']['p_' + m], 2000)
            if rule == 'A':
                add('Table 5', 'No filter, 50 label-free representatives', s, m, R[s]['noFilter_reps']['p_' + m], 2000)
                add('Table 5', 'Random forest, out-of-bag', s, m, R[s]['randomForest_oob']['p_' + m], 500)
out = []
for it in ITEMS:
    k = k_of(it['p'], it['n']); lo, hi = cp(k, it['n'])
    thr = [t for t in THR if lo <= t <= hi]
    old = lo <= THR_OLD <= hi
    if thr or old:
        out.append(dict(it, k=k, ci_low=lo, ci_high=hi, threshold=' and '.join('%.4g' % t for t in thr),
                        contains_0_0125=old, in_table_S10=bool(thr)))
json.dump(dict(n_checked=len(ITEMS), thresholds=list(THR), borderline=out), open(OUT, 'w'), indent=1, ensure_ascii=False)
n10 = sum(o['in_table_S10'] for o in out)
print('checked:', len(ITEMS), '| interval contains 0.05 or 0.0167:', n10,
      '| contains 0.0125 only (dropped):', sum(1 for o in out if not o['in_table_S10']))
for o in out:
    print('  %-3s %-8s %-48s %-14s %-4s p = %.4f  k = %3d/%d  95 %% CI %.4f-%.4f  (%s)%s' % (
        'S10' if o['in_table_S10'] else '---', o['Source'], o['Model'], o['Endpoint'], o['Metric'], o['p'], o['k'], o['n'],
        o['ci_low'], o['ci_high'], o['threshold'] or '-', '  [0.0125]' if o['contains_0_0125'] else ''))
