# -*- coding: utf-8 -*-
"""Permutation draws of the Route A1 top-15 models (Table S17: top 15 by effect size, Table S18: top 15 by VIP with one
component) from ../../code/prun_a1.py, compared with the archived objects.

For every strain the draws of a run of prun_a1.py (mode e or v1) are compared with
  - the published summaries in ../../results/revision/cache_topk_seeds/<SHORT>_A1_top15.json (e) or
    <SHORT>_A1_viptop1.json (v1): observed Q2 and AUC, p-values, median null Q2 and mean null AUC,
  - for v1 also ../../results/raw_permutations/vip1lv_final.npy,
  - the older archived draws perm_a1e_<SHORT>.npz or perm_a1vip1_<SHORT>.npz in ../../results/raw_permutations/,
    which were computed before the fallback correction (../../CHANGELOG_revision.md, section 1).
Every draw that differs from the older archived one is computed again twice, with the corrected code and with a copy of
the model function whose fallback predicts ytr.mean() as before the correction.

Usage: python a1_top15_draws_check.py <archive_root> <run_dir> <e|v1> <out.json>
<run_dir> holds perm_a1<mode>_<SHORT>.npz for HI, SA, SPneu, SPyo and PA, for example results/recalc after the commands
of ../../README.md, or ../../results/raw_permutations/a1_top15_corrected.
"""
import json, os, sys
import numpy as np
ROOT, RUN, MODE, OUT = sys.argv[1:5]
assert MODE in ('e', 'v1'), MODE
sys.path.insert(0, os.path.join(ROOT, 'code'))
import nested_generic as ng, altorder as ao
YCOL = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPneu': 'S. Pneu (S. pneumoniae)',
        'SPyo': 'S. pyo (S. pyogenes)', 'PA': 'P. A. (P. aeruginosa)'}
FUNC, NEW, OLD, SEED, OLDFILE, CACHE = {
    'e': ('alt_nested_cv', "yhat[F['i']] = 0.5", "yhat[F['i']] = ytr.mean()", 59, 'perm_a1e_%s.npz', '%s_A1_top15.json'),
    'v1': ('alt_viptop', "yhat[F['i']]=0.5", "yhat[F['i']]=ytr.mean()", 137, 'perm_a1vip1_%s.npz', '%s_A1_viptop1.json')}[MODE]
lines = open(os.path.join(ROOT, 'code', 'altorder.py'), encoding='utf-8').read().splitlines()
i0 = next(i for i, l in enumerate(lines) if l.startswith('def ' + FUNC + '('))
i1 = next((i for i in range(i0 + 1, len(lines)) if lines[i].startswith('def ') or lines[i].startswith('if __name__')), len(lines))
body = '\n'.join(lines[i0:i1])
assert body.count(NEW) == 2, body.count(NEW)
ns = dict(vars(ao)); exec(body.replace(NEW, OLD).replace('def ' + FUNC + '(', 'def old_model('), ns); old_model = ns['old_model']
new_model = getattr(ao, FUNC)
def model(f, D, groups, yv):
    return f(D, groups, yv, 1, top_k=15) if MODE == 'e' else f(D, groups, yv, 1, 15, nlv_vip=1)
def p(v, o): return float((1 + np.sum(v >= o - 1e-9)) / (1 + len(v)))
RP = os.path.join(ROOT, 'results', 'raw_permutations')
V1 = np.load(os.path.join(RP, 'vip1lv_final.npy'), allow_pickle=True).item() if MODE == 'v1' else None
res = {}
for s in YCOL:
    n = np.load(os.path.join(RUN, 'perm_a1%s_%s.npz' % (MODE, s)))
    q, a, oq, oa = n['q2s'], n['aucs'], float(n['obs_q2']), float(n['obs_auc'])
    c = json.load(open(os.path.join(ROOT, 'results', 'revision', 'cache_topk_seeds', CACHE % s)))
    run = dict(n=int(len(q)), q2=oq, auc=oa, p_q2=p(q, oq), p_auc=p(a, oa), null_q2_median=float(np.median(q)),
               null_auc_mean=float(np.mean(a)))
    r = dict(run=run, equals_published_summary=all(run[k] == c[k] for k in run if k != 'n') and run['n'] == c['nperm'])
    if V1 is not None:
        v = V1[{'HI': 'H.I.', 'SA': 'S.A.', 'SPneu': 'S.Pneu', 'SPyo': 'S.Pyo', 'PA': 'P.A.'}[s]]
        r['equals_vip1lv_final'] = bool(float(v['nq']) == oq and float(v['na']) == oa and float(v['pq']) == run['p_q2']
                                        and float(v['pa']) == run['p_auc'])
    z = np.load(os.path.join(RP, OLDFILE % s))
    zq, za = z['q2s'], z['aucs']
    idx = np.where((zq != q) | (za != a))[0]
    r['older_draws'] = dict(file=OLDFILE % s, observed_equal=bool(float(z['obs_q2']) == oq and float(z['obs_auc']) == oa),
                            p_q2=p(zq, float(z['obs_q2'])), p_auc=p(za, float(z['obs_auc'])), differing_draws=int(len(idx)))
    D = ng.build('Cleaned_Steps1-4', ycol=YCOL[s]); groups = ao.get_groups(D['X'])
    new_ok = old_ok = 0
    for k in idx:
        yv = np.random.default_rng([SEED, int(k)]).permutation(D['Y'])
        qn, an, _ = model(new_model, D, groups, yv); qo, ao_, _ = model(old_model, D, groups, yv)
        new_ok += bool(qn == q[k] and an == a[k]); old_ok += bool(qo == zq[k] and ao_ == za[k])
    r['older_draws'].update(corrected_code_reproduces_run=new_ok, former_fallback_reproduces_older=old_ok)
    res[s] = r; print(s, json.dumps(r), flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
