# -*- coding: utf-8 -*-
"""Comparison of the p-values before and after the treatment of ties (revision of 03.10.2026).

Before: result files with corrected fallback branch, strict comparison (backup in
out/strict_backup_20261003). After: the same runs, same seeds, ties counted.
Also checks that the point estimates and the Q2 p-values are unchanged (determinism)
and that every new p-value is >= the old one. Writes out/tie_effect.json and .xlsx.
"""
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)
import json, glob, os
import pandas as pd

B = (_OUT + '/strict_backup_20261003/')
if not os.path.isdir(B):   # without the backup the comparison would be incomplete and would overwrite tie_effect.json
    raise SystemExit('compare_tie.py needs ' + B + ' (results before the treatment of ties), which is not part of '
                     'this archive. Nothing was written.')
N = (_OUT + '/')
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
        fo, fn = B + 'lone/' + nm, (_OUT + '/lone/') + nm
        if os.path.exists(fo) and os.path.exists(fn):
            walk(json.load(open(fo)), json.load(open(fn)), ['lone/' + nm])
        else:
            missing.append('lone/' + nm)

import numpy as np
RP = (_ROOT + '/results/raw_permutations/')
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
print('p-values compared:', len(pv), '| changed:', summary['n_changed'], '| largest change: %.4f' % summary['max_change'])
print('moved above 0.05:', [(c['path'], round(c['old'], 4), round(c['new'], 4)) for c in summary['cross_05']])
print('moved above 0.0125:', [(c['path'], round(c['old'], 4), round(c['new'], 4)) for c in summary['cross_bonf']])
print('point estimates changed:', len(bad_est), '| Q2-p changed:', len(bad_q2), '| p became smaller:', len(neg))
print('missing:', missing)
