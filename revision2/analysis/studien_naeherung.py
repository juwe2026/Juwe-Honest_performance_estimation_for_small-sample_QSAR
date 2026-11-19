"""Post hoc, 05.11.2026: parameter-free approximation F ~ sqrt(1 + k (m-1) rbar2) for the three data sets of the study,
with the code of Revision 1 (groups altorder.build_groups, label-free representatives unsup_rep.centrality_reps),
so that prediction and observed factor (Additional file 2) come from the same implementation.
Usage: QSAR_ROUTEA_XLSX=<RouteA_matrix.xlsx> QSAR_RECALC_DIR=<tmp> python -I studien_naeherung.py <archive_root> <out.json>"""
import sys, os, json, numpy as np, pandas as pd
from scipy.stats import norm, rankdata
AF7 = sys.argv[1]; sys.path.insert(0, os.path.join(AF7, 'code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur, gaucher_pipe as gp
q = 0.05; c = norm.ppf(1 - q / 2); k = 2 * c ** 2 * norm.pdf(c) ** 2 / (q * (1 - q))
def rbar2(X, idx=None):
    R = rankdata(X, axis=0) if idx is None else rankdata(X[:, idx], axis=0)
    C = np.nan_to_num(np.corrcoef(R, rowvar=False)); iu = np.triu_indices(C.shape[0], 1)
    return float(np.mean(C[iu] ** 2))
def lauf(name, X, G, reps):
    p = X.shape[1]; a = rbar2(X); r = rbar2(X, np.asarray(reps))
    return dict(p=p, G=len(G), rbar2_alle=a, rbar2_reps=r,
                F_A2_vorhergesagt=float(np.sqrt(1 + k * (p - 1) * a)), F_A1u_vorhergesagt=float(np.sqrt(1 + k * (len(G) - 1) * r)))
out = {'k': k}
D = ng.build('Cleaned_Steps1-4', ycol='H. I. (H. influenzae)'); Xm = D['X']
Gm, _ = ao.build_groups(Xm); out['Monoterpenoide'] = lauf('Mono', Xm, Gm, ur.centrality_reps(Xm, Gm))
Xp = pd.read_csv(os.path.join(AF7, 'revision/lone/desc_pool.csv')); Xp = Xp[[c_ for c_ in Xp.columns if c_ != 'No']].values.astype(float)
Gp, _ = ao.build_groups(Xp); out['Octapeptide'] = lauf('Pep', Xp, Gp, ur.centrality_reps(Xp, Gp))
X, y, mz, _ = gp.load(os.path.join(AF7, 'data/Gaucherdata_paper39_median_normalised.csv')); Xg, _, _ = gp.clean(X, mz)
Gg = gp.build_groups(Xg); Gg = Gg[0] if isinstance(Gg, tuple) else Gg
out['Gaucher'] = lauf('Gaucher', Xg, Gg, gp.centrality_reps(Xg, Gg))
json.dump(out, open(sys.argv[2], 'w'), indent=1); print(json.dumps(out, indent=1))
