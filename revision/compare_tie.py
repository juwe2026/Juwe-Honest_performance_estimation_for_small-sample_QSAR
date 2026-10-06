# -*- coding: utf-8 -*-
"""Vergleich der p-Werte vor und nach der Gleichstandsbehandlung (Revision 03.10.2026).

Vorher: Ergebnisdateien mit korrigiertem Rueckfallzweig, strikter Vergleich (Sicherung in
out/strict_backup_20261003). Nachher: dieselben Laeufe, gleiche Seeds, Gleichstaende mitgezaehlt.
Prueft ausserdem, dass die Punktschaetzer und die Q2-p-Werte unveraendert sind (Determinismus)
und dass jeder neue p-Wert >= dem alten ist. Schreibt out/tie_effect.json und .xlsx.
"""
import json, glob, os
import pandas as pd

B = '/home/claude/rev/out/strict_backup_20261003/'
N = '/home/claude/rev/out/'
rows = []

def walk(old, new, path):
    if isinstance(old, dict) and isinstance(new, dict):
        keys = set(old) & set(new)
        for k in sorted(keys, key=str):
            walk(old[k], new[k], path + [str(k)])
        return
    if isinstance(old, (int, float)) and isinstance(new, (int, float)):
        rows.append(dict(path='/'.join(path), old=float(old), new=float(new)))

PAIRS = [('I_fallback_seeds.json', 'I_fallback_seeds.json'), ('A2_HI_seed42.json', 'A2_HI_seed42.json'),
         ('A_nofilter.json', 'A_nofilter.json'), ('D_rf.json', 'D_rf.json'), ('D_rf_morgan.json', 'D_rf_morgan.json'),
         ('E_misc.json', 'E_misc.json'), ('morgan2_fixed.json', 'morgan2_fixed.json'), ('C_routeB.json', 'C_routeB.json'),
         ('H_median_fixed_HI.json', 'H_median_fixed_HI.json'), ('F_no3d_fixed.json', 'F_no3d_fixed.json')]
for fo in sorted(glob.glob(B + 'cache_topk_seeds/*.json')):
    PAIRS.append(('cache_topk_seeds/' + os.path.basename(fo), 'cache_topk_seeds/' + os.path.basename(fo)))
missing = []
for a, b in PAIRS:
    if not (os.path.exists(B + a) and os.path.exists(N + b)):
        missing.append(a); continue
    walk(json.load(open(B + a)), json.load(open(N + b)), [a])
for r_ in ('A', 'C'):
    for nm, sub in (('res_%s.json' % r_, ''), ('res_%s_fix.json' % r_, '')):
        fo, fn = B + 'lone/' + nm, '/home/claude/rev/lone/' + nm
        if os.path.exists(fo) and os.path.exists(fn):
            walk(json.load(open(fo)), json.load(open(fn)), ['lone/' + nm])
        else:
            missing.append('lone/' + nm)

import numpy as np
RP = '/home/claude/rev/QSAR_selection_leakage/results/raw_permutations/'
for pre in ('permB_12', 'permB_12nlv3', 'permB_25', 'gperm_A1', 'gperm_A1u', 'gperm_A2', 'gperm2_A1', 'gperm2_A1u', 'gperm2_A2'):
    for s_ in (('HI', 'SA', 'SPneu', 'SPyo') if pre.startswith('permB') else ('',)):
        z = np.load(RP + (pre + '_' + s_ if s_ else pre) + '.npz')
        for m, nk, ok in (('p_q2', 'q2s', 'obs_q2'), ('p_auc', 'aucs', 'obs_auc')):
            nul, ob = z[nk], float(z[ok])
            rows.append(dict(path='archive/%s/%s/%s' % (pre, s_ or '-', m), old=float((1 + np.sum(nul >= ob)) / (1 + len(nul))),
                             new=float((1 + np.sum(nul >= ob - 1e-9)) / (1 + len(nul)))))
df = pd.DataFrame(rows)
leaf = df.path.str.split('/').str[-1]
isp = leaf.isin(['p_q2', 'p_auc', 'p'])
ign = df.path.str.contains('seconds|_old|hits/', regex=True)
pv = df[isp & ~ign].copy()
pv['metric'] = leaf[isp & ~ign]
pv['diff'] = pv.new - pv.old
est = df[~isp & ~ign & ~leaf.isin(['nperm', 'blowups'])]
bad_est = est[(est.new - est.old).abs() > 1e-9]
bad_q2 = pv[(pv.metric == 'p_q2') & (pv['diff'].abs() > 1e-12)]
neg = pv[pv['diff'] < -1e-12]
cross = pv[(pv.old < 0.05) & (pv.new >= 0.05)]
cross_b = pv[(pv.old < 0.0125) & (pv.new >= 0.0125)]
summary = dict(n_p=len(pv), n_changed=int((pv['diff'].abs() > 1e-12).sum()), max_change=float(pv['diff'].max()),
               changed=pv[pv['diff'].abs() > 1e-12][['path', 'old', 'new']].to_dict('records'),
               cross_05=cross[['path', 'old', 'new']].to_dict('records'),
               cross_bonf=cross_b[['path', 'old', 'new']].to_dict('records'),
               estimate_changes=bad_est[['path', 'old', 'new']].to_dict('records')[:50],
               q2_p_changes=bad_q2[['path', 'old', 'new']].to_dict('records'),
               negative=neg[['path', 'old', 'new']].to_dict('records'), missing=missing)
json.dump(summary, open(N + 'tie_effect.json', 'w'), indent=1)
pv.to_excel(N + 'tie_effect.xlsx', index=False)
print('p-Werte verglichen:', len(pv), '| geaendert:', summary['n_changed'], '| groesste Aenderung: %.4f' % summary['max_change'])
print('ueber 0,05 gerutscht:', [(c['path'], round(c['old'], 4), round(c['new'], 4)) for c in summary['cross_05']])
print('ueber 0,0125 gerutscht:', [(c['path'], round(c['old'], 4), round(c['new'], 4)) for c in summary['cross_bonf']])
print('Punktschaetzer veraendert:', len(bad_est), '| Q2-p veraendert:', len(bad_q2), '| p kleiner geworden:', len(neg))
print('fehlend:', missing)
