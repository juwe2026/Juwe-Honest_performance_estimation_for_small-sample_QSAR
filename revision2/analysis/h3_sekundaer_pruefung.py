"""Post hoc descriptive check (Revision 2):
Why does the calibrated inflation factor correlate POSITIVELY with G/p under A1 and A1u,
but negatively under A2 (H3)?

Decomposition: for m tests with hit rate q and mean pairwise correlation rho of the
hit indicators, Var(hit count) = m q (1-q) [1 + (m-1) rho], hence
F^2 = 1 + (m-1) rho  (F = calibrated inflation factor).
Under A1/A1u m = G (number of groups), under A2 m = p (all descriptors).
The check asks whether the number of tests (m) or the dependence per pair of tests (rho) carries
the association with G/p. For this, directly from the descriptor matrices: mean r^2
between the label-free representatives (A1u) and between all descriptors (A2).
Usage: python -I h3_sekundaer_pruefung.py <datasets_k> <results.xlsx> <scripts> <out.json>
"""
import sys, json, numpy as np, pandas as pd, openpyxl
from scipy.stats import spearmanr
sys.path.insert(0, sys.argv[3])
from routen import zentralitaets_reps

DSDIR, XLSX, OUT = sys.argv[1], sys.argv[2], sys.argv[4]
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
rows = list(wb['Summary'].iter_rows(values_only=True))
h = [i for i, r in enumerate(rows) if r and r[0] == 'Data set'][0]; H = list(rows[h])
S = {r[0]: dict(zip(H, r)) for r in rows[h+1:] if r and r[0] and str(r[0]).startswith('DS') and int(r[0][2:]) <= 16}

def mittel_r2(A, idx):
    C = np.corrcoef(A[:, idx], rowvar=False); iu = np.triu_indices(len(idx), 1)
    return float(np.mean(np.nan_to_num(C[iu]) ** 2)), float(np.mean(np.abs(np.nan_to_num(C[iu]))))

res = []
for ds in sorted(S):
    Xdf = pd.read_csv(f'{DSDIR}/{ds}_X.csv'); cols = list(Xdf.columns); A = Xdf.to_numpy(float)
    info = json.load(open(f'{DSDIR}/{ds}_info.json'))
    gr = [[cols.index(c) for c in g if c in cols] for g in info['gruppen']]   # stored groups of the main analysis
    reps = zentralitaets_reps(A, gr)                                           # as run_benchmark_k.lade()
    r2_reps, abs_reps = mittel_r2(A, reps)
    r2_alle, abs_alle = mittel_r2(A, np.arange(A.shape[1]))
    s = S[ds]; G = s['G']; p = s['p (cleaned)']
    d = dict(ds=ds, n=s['n'], p=p, G=G, G_neu=len(gr), Gp=s['G/p'],
             F_A1u=s['Calibrated inflation A1u'], F_A1=s['Calibrated inflation A1'], F_A2=s['Calibrated inflation A2'],
             r2_reps=r2_reps, absr_reps=abs_reps, r2_alle=r2_alle, absr_alle=abs_alle)
    for k, m in [('A1u', G), ('A1', G), ('A2', p)]:
        d[f'rho_{k}'] = (d[f'F_{k}'] ** 2 - 1) / (m - 1)       # implied hit correlation per pair
    res.append(d)

df = pd.DataFrame(res)
print(df[['ds','n','p','G','G_neu','Gp','F_A1u','F_A1','F_A2','r2_reps','r2_alle','rho_A1u','rho_A1','rho_A2']].round(4).to_string(index=False))
def sp(a, b):
    r = spearmanr(df[a], df[b]); return dict(rho=round(float(r.statistic), 3), p_zweiseitig=round(float(r.pvalue), 4))
tests = {}
for a, b in [('F_A2','Gp'), ('F_A1','Gp'), ('F_A1u','Gp'),
             ('F_A1','G'), ('F_A1u','G'), ('F_A2','p'), ('G','Gp'), ('p','Gp'),
             ('rho_A2','Gp'), ('rho_A1','Gp'), ('rho_A1u','Gp'),
             ('r2_reps','Gp'), ('r2_alle','Gp'), ('rho_A1u','r2_reps'), ('rho_A2','r2_alle')]:
    tests[f'{a} ~ {b}'] = sp(a, b)
for k, v in tests.items(): print(f'{k:22s} rho {v["rho"]:+.3f}  p {v["p_zweiseitig"]:.4f}')
# How well does m alone explain the factor? F predicted with constant rho (median)
for k, m in [('A1u','G'), ('A1','G'), ('A2','p')]:
    rho_med = float(np.median(df[f'rho_{k}']))
    pred = np.sqrt(1 + (df[m] - 1) * rho_med)
    r = spearmanr(pred, df[f'F_{k}']).statistic
    tests[f'F_{k} vorhergesagt aus m={m} bei konstantem rho'] = dict(rho_median=round(rho_med, 5), spearman_mit_beobachtet=round(float(r), 3))
    print(k, 'rho_median', round(rho_med, 5), 'Spearman(predicted, observed)', round(float(r), 3))
json.dump(dict(tabelle=df.round(6).to_dict(orient='records'), tests=tests), open(OUT, 'w'), indent=1)
