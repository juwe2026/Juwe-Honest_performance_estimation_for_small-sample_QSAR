# -*- coding: utf-8 -*-
"""Table S5 for Revision 2: the 100 random tie-break rules of the Route A1u grouping, evaluated against 0.0167.

The counting step of jobB_hits.py is repeated with the same seeds (null permutations default_rng([8202, j]),
j = 0..999, random tie-breaks default_rng(8203), 100 groupings), the same filter (|r_rb| > 0.3 and Mann-Whitney
p <= 0.05) and the same p-value (pval_tie.py, ties within 1e-9 counted). New is only that the p-value of every
single rule is kept, so that the share of rules with p <= 0.0167 can be reported. Before writing, the script checks
that the published rule, the first-index rule and the summary of the 100 random rules reproduce B_hits.json exactly.

Usage: python revision2/analysis/jobB_tiebreak_0167.py ARCHIVE_ROOT OUT.json   (needs base.pkl, see make_base.py and common.py)
"""
import sys, os, json
ROOT, OUTF = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(ROOT, 'revision'))
from common import *                      # base.pkl, X, ADJ, Y, KEYS, N, pval (pval_tie.py)
import nested_generic as ng
from scipy.stats import rankdata

R = rankdata(X, axis=0); p = X.shape[1]
Cp = np.abs(np.nan_to_num(np.corrcoef(X, rowvar=False))); Cs = np.abs(np.nan_to_num(np.corrcoef(R, rowvar=False)))
C = np.maximum(Cp, Cs); np.fill_diagonal(C, np.nan)

def groups(tb, rng=None):                 # identical to jobB_hits.py
    rem = np.ones(p, bool); deg = ADJ.sum(1).astype(int); Gs = []; var = X.var(0)
    while rem.any():
        d = np.where(rem, deg, -1); mx = d.max()
        if mx <= 0: Gs += [[int(i)] for i in np.where(rem)[0]]; break
        cand = np.where(d == mx)[0]
        hub = cand[np.argmax(var[cand])] if tb == 'var' else (cand[0] if tb == 'first' else rng.choice(cand))
        mem = np.unique(np.concatenate([[hub], np.where(ADJ[hub] & rem)[0]])); Gs.append(mem.tolist())
        for m_ in mem: deg -= ADJ[m_].astype(int)
        rem[mem] = False
    return Gs
def creps(Gs):
    return np.array([g[0] if len(g) == 1 else g[int(np.nanargmax(np.nanmean(C[np.ix_(g, g)], 1)))] for g in Gs])

GV = {'variance (as published)': groups('var'), 'first index': groups('first')}
rng = np.random.default_rng(8203)
for r in range(100): GV['random_%03d' % r] = groups('rand', rng)
RP = {g: creps(v) for g, v in GV.items()}
assert len(GV['variance (as published)']) == 182
NP = 1000
pass_obs = {}; nulls = {g: {k: [] for k in KEYS} for g in GV}
for k in KEYS:
    y = Y[k]; eff, pv = ng.eff_filter(R, y, N); pass_obs[k] = (eff > 0.3) & (pv <= 0.05)
    for j in range(NP):
        yp = np.random.default_rng([8202, j]).permutation(y); e, pv_ = ng.eff_filter(R, yp, N); ps = (e > 0.3) & (pv_ <= 0.05)
        for g in GV: nulls[g][k].append(int(ps[RP[g]].sum()))
res = {}
for g in GV:
    res[g] = {'n_groups': len(GV[g])}
    for k in KEYS:
        o = int(pass_obs[k][RP[g]].sum()); nl = np.array(nulls[g][k])
        res[g][k] = dict(obs=o, null_mean=float(nl.mean()), null_sd=float(nl.std(ddof=1)), p=pval(nl, o))
rand = [g for g in GV if g.startswith('random')]

# reproduction check against the archived result
B = json.load(open(os.path.join(ROOT, 'results', 'revision', 'B_hits.json')))
for g in ('variance (as published)', 'first index'):
    for k in KEYS:
        a, b = B['tiebreak_fixed'][g][k], res[g][k]
        assert a['obs'] == b['obs'] and abs(a['p'] - b['p']) < 1e-12, (g, k, a, b)
summ = {}
for k in KEYS:
    ps = np.array([res[g][k]['p'] for g in rand]); ob = np.array([res[g][k]['obs'] for g in rand])
    a = B['random_summary'][k]
    assert (a['obs_min'], a['obs_max']) == (int(ob.min()), int(ob.max())) and abs(a['p_median'] - float(np.median(ps))) < 1e-12
    assert abs(a['frac_p_le_0125'] - float(np.mean(ps <= 0.0125))) < 1e-12 and abs(a['frac_p_le_05'] - float(np.mean(ps <= 0.05))) < 1e-12
    summ[k] = dict(obs_min=int(ob.min()), obs_max=int(ob.max()), obs_median=float(np.median(ob)), p_min=float(ps.min()),
                   p_max=float(ps.max()), p_median=float(np.median(ps)), frac_p_le_05=float(np.mean(ps <= 0.05)),
                   frac_p_le_0167=float(np.mean(ps <= 0.05 / 3)), frac_p_le_0125=float(np.mean(ps <= 0.0125)))
assert [res[g]['n_groups'] for g in rand] == B['random_ngroups']
json.dump(dict(reproduces_B_hits=True, summary=summ, per_rule={g: {k: res[g][k] for k in KEYS} | {'n_groups': res[g]['n_groups']} for g in rand},
               seeds=B['seeds']), open(OUTF, 'w'), indent=1)
for k in KEYS:
    s = summ[k]
    print('%-6s p median %.4f (%.4f-%.4f) | share p<=0.05 %3.0f %% | p<=0.0167 %3.0f %% | p<=0.0125 %3.0f %%' % (
        k, s['p_median'], s['p_min'], s['p_max'], 100 * s['frac_p_le_05'], 100 * s['frac_p_le_0167'], 100 * s['frac_p_le_0125']))
print('reproduces B_hits.json: yes')
