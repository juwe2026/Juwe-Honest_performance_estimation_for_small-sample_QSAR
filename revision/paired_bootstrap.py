# Paired bootstrap of differences in honest Q2 between models (revision, October 2026): Table S22 and the Methods
# paragraph on model comparisons.
#
# Input: the observed leave-one-out (or out-of-bag) predictions of every compared model, written by
# bootstrap_predictions_ds1.py (monoterpenoids) and bootstrap_predictions_ds2.py (octapeptides, main activity definition).
# The predictions are held fixed, nothing is refitted.
#
# For every endpoint, 20,000 bootstrap resamples of the compound-prediction pairs are drawn, stratified by class
# (actives and non-actives resampled separately with replacement), with numpy.random.default_rng(8302). The same
# resamples are used for every model of that endpoint, so the two Q2 values of a comparison are paired. The first 2,000
# resamples are those of Table 1, panel d (jobE_misc.py). For each comparison the script reports Q2 of model 1 and of
# model 2, their difference, its 95 % percentile interval, and a two-sided bootstrap p-value
#   p = 2 * min((#{d <= 0} + 1) / (B + 1), (#{d >= 0} + 1) / (B + 1)),  capped at 1,
# and applies the Holm and the Benjamini-Hochberg procedure (5 % each) over all 104 comparisons.
#
# Usage: python revision/paired_bootstrap.py   (paths as in the other revision scripts, see README_revision.md)
import json, csv
import numpy as np

RES = '/home/claude/rev/out/'          # results/revision/ in this archive
B = 20000
SEED = 8302

def q2(y, yh):
    return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)

def resamples(y):
    rng = np.random.default_rng(SEED); i1 = np.where(y == 1)[0]; i0 = np.where(y == 0)[0]
    return [np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]) for _ in range(B)]

# model keys of the prediction files and the labels printed in Table S22
LAB = {'NoSel2208': 'PLS-DA, all 2,208 descriptors', 'RF2208': 'Random forest, all 2,208', 'A1': 'Route A1', 'A2': 'Route A2',
       'A1u': 'Route A1u', 'B12': 'Route B (12, 1 LV)', 'A1_top15': 'Route A1 top-15', 'A2_top15': 'Route A2 top-15',
       'A1_top6': 'Route A1 top-6', 'A1_vip15': 'Route A1 top-15, VIP', 'A2_vip15': 'Route A2 top-15, VIP',
       'B25': '25 literature descriptors', 'B11_noPhenol': 'Route B without Phenol_flag (11)', 'B12_3LV': 'Route B, 3 LV',
       'NoSel1275': 'PLS-DA, all 1,275 descriptors', 'RF1275': 'Random forest, all 1,275', 'Reps50': '50 representatives, no filter'}
# (model 1, model 2, label of model 1 if different, label of model 2 if different); difference = model 1 minus model 2
DS1 = [('NoSel2208', 'A1'), ('NoSel2208', 'A2'), ('RF2208', 'A1'), ('RF2208', 'A2'), ('B12', 'A1'),
       ('B12', 'NoSel2208', None, 'PLS-DA, all 2,208'), ('B12', 'RF2208'), ('A1', 'A2'), ('A1', 'A1u'), ('A2', 'A1u'),
       ('A1_top15', 'A1'), ('A2_top15', 'A2'), ('A1_top6', 'A1'),
       ('A1_top15', 'A1_vip15', 'Route A1 top-15, effect size'), ('A2_top15', 'A2_vip15', 'Route A2 top-15, effect size'),
       ('B12', 'B25', 'Route B, 12 descriptors'), ('B12', 'B11_noPhenol', 'Route B, 12 descriptors'),
       ('B12_3LV', 'B12', None, 'Route B, 1 LV')]
DS2 = [('NoSel1275', 'A1'), ('NoSel1275', 'A2'), ('NoSel1275', 'A1u'), ('RF1275', 'A1'), ('RF1275', 'A1u'), ('Reps50', 'A1'),
       ('A1', 'A2'), ('A1', 'A1u')]
EP = {'ds1': [('HI', 'H. influenzae'), ('SA', 'S. aureus'), ('SPneu', 'S. pneumoniae'), ('SPyo', 'S. pyogenes')],
      'ds2': [('EC', 'E. coli'), ('SA', 'S. aureus'), ('PA', 'P. aeruginosa'), ('AB', 'A. baumannii')]}
PRED = {'ds1': json.load(open(RES + 'bootstrap_predictions_ds1.json')), 'ds2': json.load(open(RES + 'bootstrap_predictions_ds2.json'))}

rows = []
for ds, comps in (('ds1', DS1), ('ds2', DS2)):
    for k, ep in EP[ds]:
        y = np.array(PRED[ds][k]['y']); P = {m: np.array(v) for m, v in PRED[ds][k]['pred'].items()}
        I = resamples(y)
        qb = {}
        for c in comps:
            for m in c[:2]:
                if m not in qb: qb[m] = np.array([q2(y[i], P[m][i]) for i in I])
        for c in comps:
            m1, m2 = c[:2]
            l1 = c[2] if len(c) > 2 and c[2] else LAB[m1]; l2 = c[3] if len(c) > 3 and c[3] else LAB[m2]
            d = qb[m1] - qb[m2]; d0 = q2(y, P[m1]) - q2(y, P[m2])
            lo, hi = np.percentile(d, [2.5, 97.5])
            p = min(1.0, 2 * min((np.sum(d <= 0) + 1) / (B + 1), (np.sum(d >= 0) + 1) / (B + 1)))
            rows.append(dict(data_set=ds, model_1=l1, model_2=l2, endpoint=ep, q2_model_1=q2(y, P[m1]), q2_model_2=q2(y, P[m2]),
                             difference=d0, ci95_lower=lo, ci95_upper=hi, ci_excludes_zero=bool(lo > 0 or hi < 0), p_two_sided=p))
# sort rows as printed: comparison by comparison, endpoints in the order above
order = {(ds, i): n for n, (ds, i) in enumerate([('ds1', i) for i in range(len(DS1))] + [('ds2', i) for i in range(len(DS2))])}
def key(r):
    comps = DS1 if r['data_set'] == 'ds1' else DS2
    ci = [n for n, c in enumerate(comps) if (c[2] if len(c) > 2 and c[2] else LAB[c[0]]) == r['model_1']
          and (c[3] if len(c) > 3 and c[3] else LAB[c[1]]) == r['model_2']][0]
    return (order[(r['data_set'], ci)], [e for _, e in EP[r['data_set']]].index(r['endpoint']))
rows.sort(key=key)
assert len(rows) == 104
# Holm and Benjamini-Hochberg over all 104 comparisons, 5 % each
p = np.array([r['p_two_sided'] for r in rows]); m = len(p); o = np.argsort(p, kind='mergesort')
holm = np.zeros(m, bool)
for j, i in enumerate(o):
    if p[i] <= 0.05 / (m - j): holm[i] = True
    else: break
ok = [j for j, i in enumerate(o) if p[i] <= 0.05 * (j + 1) / m]
bh = np.zeros(m, bool)
if ok: bh[o[:max(ok) + 1]] = True
for r, h, b in zip(rows, holm, bh): r['holm_5pct'] = bool(h); r['benjamini_hochberg_5pct'] = bool(b)
n_sig = sum(r['ci_excludes_zero'] for r in rows)
print('%d comparisons, %d intervals exclude zero (about %.1f expected by chance), Holm %d, Benjamini-Hochberg %d'
      % (m, n_sig, 0.05 * m, holm.sum(), bh.sum()))
json.dump(dict(B=B, seed=SEED, n_comparisons=m, n_ci_excludes_zero=n_sig, n_holm=int(holm.sum()), n_bh=int(bh.sum()), rows=rows),
          open(RES + 'paired_bootstrap_Q2.json', 'w'), indent=1, default=float)
with open(RES + 'paired_bootstrap_Q2.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
