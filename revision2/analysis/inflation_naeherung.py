"""Post hoc, 05.11.2026: parameter-free approximation of the calibrated inflation factor of the hit count.
F ~ sqrt(1 + k (m-1) rbar2), k = 2 c^2 phi(c)^2 / (q (1-q)), c = Phi^-1(1 - q/2), q = 0.05 (nominal, not fitted).
Leading term of the Hermite expansion for the correlation of two two-sided hit indicators. rbar2 = mean square
of the rank correlation (Spearman) between the tested descriptors, because the Mann-Whitney test is rank-based.
For comparison also with Pearson. m = p (Route A2) or G (A1u, A1, representatives chosen label-free by centrality).
Usage: python -I inflation_naeherung.py <h3_ergebnis.json> <datasets_k> <benchmark_scripts> <out.json>"""
import sys, json, numpy as np, pandas as pd
from scipy.stats import spearmanr, norm, rankdata
sys.path.insert(0, sys.argv[3]); from routen import zentralitaets_reps
d = pd.DataFrame(json.load(open(sys.argv[1]))['tabelle']); D = sys.argv[2]
rs_all, rs_rep = [], []
for ds in d.ds:
    Xdf = pd.read_csv(f'{D}/{ds}_X.csv'); cols = list(Xdf.columns); A = Xdf.to_numpy(float)
    info = json.load(open(f'{D}/{ds}_info.json')); gr = [[cols.index(c) for c in g if c in cols] for g in info['gruppen']]
    reps = zentralitaets_reps(A, gr)
    C = np.nan_to_num(np.corrcoef(rankdata(A, axis=0), rowvar=False)); iu = np.triu_indices(A.shape[1], 1)
    rs_all.append(float(np.mean(C[iu] ** 2))); Cr = C[np.ix_(reps, reps)]; ir = np.triu_indices(len(reps), 1); rs_rep.append(float(np.mean(Cr[ir] ** 2)))
d['rs2_alle'] = rs_all; d['rs2_reps'] = rs_rep
q = 0.05; c = norm.ppf(1 - q / 2); k = 2 * c ** 2 * norm.pdf(c) ** 2 / (q * (1 - q))
out = dict(k=k, serien={}, zusammenfassung={})
for lab, a, r in [('spearman', d.rs2_alle, d.rs2_reps), ('pearson', d.r2_alle, d.r2_reps)]:
    pr = {'A2': np.sqrt(1 + k * (d.p - 1) * a), 'A1u': np.sqrt(1 + k * (d.G - 1) * r), 'A1': np.sqrt(1 + k * (d.G - 1) * r)}
    for rt in pr:
        ratio = d[f'F_{rt}'] / pr[rt]
        out['zusammenfassung'][f'{lab}_{rt}'] = dict(spearman=round(float(spearmanr(pr[rt], d[f'F_{rt}']).statistic), 3),
            verhaeltnis_median=round(float(np.median(ratio)), 3), verhaeltnis_min=round(float(ratio.min()), 3), verhaeltnis_max=round(float(ratio.max()), 3))
        for i, ds in enumerate(d.ds): out['serien'].setdefault(ds, {})[f'{lab}_{rt}_vorhergesagt'] = round(float(pr[rt][i]), 3)
for i, ds in enumerate(d.ds):
    out['serien'][ds].update(beobachtet_A2=d.F_A2[i], beobachtet_A1u=d.F_A1u[i], beobachtet_A1=d.F_A1[i], rs2_alle=round(rs_all[i], 5), rs2_reps=round(rs_rep[i], 5))
json.dump(out, open(sys.argv[4], 'w'), indent=1); print(json.dumps(out['zusammenfassung'], indent=1))
