# -*- coding: utf-8 -*-
"""Clopper-Pearson-Intervalle (95 %) fuer die Monte-Carlo-Unsicherheit der berichteten Permutations-p-Werte.

k = Zahl der permutierten Werte >= beobachtet (mit Gleichstandsbehandlung), n = Zahl der Permutationen;
berichtet wird p = (1 + k)/(1 + n), das Intervall gilt fuer den Anteil k/n, den das Monte-Carlo-Verfahren
schaetzt. Ausgewiesen werden alle in Tables 1 bis 4 und 6 berichteten p-Werte der auswertbaren Endpunkte (Tables 1 bis 5),
deren Intervall 0,05 oder 0,0125 einschliesst (Zusage im Antwortschreiben, G2-7; Entscheidung EW 03.10.2026).
Schreibt out/cp_borderline.json.
"""
import json
import numpy as np
from scipy.stats import beta

O = '/home/claude/rev/out/'
RP = '/home/claude/rev/QSAR_selection_leakage/results/raw_permutations/'
J = lambda f: json.load(open(f))

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
      'EC': 'E. coli', 'AB': 'A. baumannii'}
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
        _b = J(O + 'B25_cineole0_seed701.json')[s]; add('Table 2', '25 descriptors', s, m, _b['p_' + m], 2000)   # G2-20: Cineol = 0
        add('Table 3', 'Route A1 top-15 by effect size, nested', s, m, J(O + 'cache_topk_seeds/%s_A1_top15.json' % s)['p_' + m], 2000)
        add('Table 4', 'Route A1 top-15 by VIP, nested', s, m, J(O + 'cache_topk_seeds/%s_A1_viptop1.json' % s)['p_' + m], 2000)
for rule, lab in (('A', ''), ('C', ', sensitivity definition')):
    R = J('/home/claude/rev/lone/res_%s.json' % rule); RX = J('/home/claude/rev/lone/res_%s_fix.json' % rule)
    for s in ('EC', 'SA', 'AB'):
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
    thr = [t for t in (0.05, 0.0125) if lo <= t <= hi]
    if thr:
        out.append(dict(it, k=k, ci_low=lo, ci_high=hi, threshold=' and '.join('%g' % t for t in thr)))
json.dump(dict(n_checked=len(ITEMS), borderline=out), open(O + 'cp_borderline.json', 'w'), indent=1)
print('geprueft:', len(ITEMS), '| Intervall schliesst 0,05 oder 0,0125 ein:', len(out))
for o in out:
    print('  %-8s %-48s %-14s %-4s p = %.4f  k = %d/%d  95 %% CI %.4f–%.4f  (%s)' % (
        o['Source'], o['Model'], o['Endpoint'], o['Metric'], o['p'], o['k'], o['n'], o['ci_low'], o['ci_high'], o['threshold']))
