"""Post hoc, 05.11.2026: check of the parameter-free approximation of the calibrated inflation factor with a finite sample.
Under label permutation the centred label vector lies in an (n-1)-dimensional space. Model: Mann-Whitney Z = sqrt(d) a'u,
u uniformly distributed on the unit sphere in d = n-1 dimensions, a = normalised centred rank vector of the descriptor.
Variants of the indicator correlation g(rho) of two descriptors with rank correlation rho (threshold c = 1.96, nothing fitted):
 V0  k rho^2                         (leading term, infinite n, previous approximation)
 V1  k (d rho^2 - 1)/(d + 2)         (leading term, finite n)
 V2  exact in the sphere model       (all orders, finite n, numerical integration)
 V3  exact bivariate normal          (all orders, infinite n)
F^2 = 1 + (m-1) * mean_{pairs} g(rho).
Usage: QSAR_ROUTEA_XLSX=<RouteA_matrix.xlsx> QSAR_RECALC_DIR=<recalc> python -I naeherung_endlich_n.py <archive_root> <h3_ergebnis.json> <datasets_k> <benchmark_scripts> <out.json>"""
import sys, os, json, numpy as np, pandas as pd
from scipy.stats import norm, rankdata, multivariate_normal, spearmanr
from scipy.special import betainc, beta as Bfn
from scipy.integrate import quad
q = 0.05; c = norm.ppf(1 - q / 2); k = 2 * c ** 2 * norm.pdf(c) ** 2 / (q * (1 - q))
RG = np.linspace(0, 0.999, 1000)

def cdf_t(x, m):  # first coordinate of u ~ uniform on S^{m-1} in R^m
    x = np.clip(x, -1, 1); return 0.5 + 0.5 * np.sign(x) * betainc(0.5, (m - 1) / 2, x ** 2)
def g_kugel(rho, d):
    h = c / np.sqrt(d); qd = 1 - betainc(0.5, (d - 1) / 2, h ** 2)
    f = lambda s: (1 - s * s) ** ((d - 3) / 2) / Bfn(0.5, (d - 1) / 2)
    def inner(s):
        w = np.sqrt(1 - rho ** 2) * np.sqrt(1 - s * s)
        return f(s) * ((1 - cdf_t((h - rho * s) / w, d - 1)) + cdf_t((-h - rho * s) / w, d - 1))
    p11 = 2 * quad(inner, h, 1, limit=200, epsabs=1e-13)[0]
    return (p11 - qd ** 2) / (qd * (1 - qd)), qd
def g_normal(rho):
    if rho == 0: return 0.0
    S = [[1, rho], [rho, 1]]; P = lambda a, b: multivariate_normal.cdf([a, b], mean=[0, 0], cov=S)
    inside = P(c, c) - P(-c, c) - P(c, -c) + P(-c, -c)  # P(|X|<c,|Y|<c)
    p11 = 1 - 2 * (1 - q) + inside
    return (p11 - q * q) / (q * (1 - q))
_cache = {}
def tabelle_kugel(d):
    if d not in _cache: _cache[d] = np.array([g_kugel(r, d)[0] for r in RG])
    return _cache[d]
TN = np.array([g_normal(r) for r in RG])
def vorhersagen(rho, d, m):
    a = np.abs(rho); r2 = float(np.mean(a ** 2))
    g = {'V0': k * r2, 'V1': k * (d * r2 - 1) / (d + 2),
         'V2': float(np.mean(np.interp(a, RG, tabelle_kugel(d)))), 'V3': float(np.mean(np.interp(a, RG, TN)))}
    return {v: float(np.sqrt(max(1 + (m - 1) * gv, 0))) for v, gv in g.items()} | {'rbar2': r2, 'mittel_g': g}
def paare(X, idx=None):
    A = X if idx is None else X[:, idx]
    C = np.nan_to_num(np.corrcoef(rankdata(A, axis=0), rowvar=False)); iu = np.triu_indices(C.shape[0], 1); return C[iu]

out = dict(k=k, c=float(c), modell='Kugelmodell d = n-1', benchmark={}, studien={})
# check of the sphere model: marginal hit rate and g(0)
out['kontrolle'] = {str(d): dict(q_d=g_kugel(0.0, d)[1], g0=g_kugel(0.0, d)[0], V1_g0=-k / (d + 2)) for d in (29, 38, 49)}
# Benchmark
sys.path.insert(0, sys.argv[4]); from routen import zentralitaets_reps
T = pd.DataFrame(json.load(open(sys.argv[2]))['tabelle']); D = sys.argv[3]
for _, r in T.iterrows():
    Xdf = pd.read_csv(f'{D}/{r.ds}_X.csv'); cols = list(Xdf.columns); A = Xdf.to_numpy(float)
    info = json.load(open(f'{D}/{r.ds}_info.json')); gr = [[cols.index(c_) for c_ in g if c_ in cols] for g in info['gruppen']]
    reps = zentralitaets_reps(A, gr); n = A.shape[0]; d = n - 1
    out['benchmark'][r.ds] = dict(n=n, p=A.shape[1], G=len(reps), beob_A2=float(r.F_A2), beob_A1u=float(r.F_A1u),
        A2=vorhersagen(paare(A), d, A.shape[1]), A1u=vorhersagen(paare(A, reps), d, len(reps)))
    print(r.ds, n, round(r.F_A2, 2), {v: round(out['benchmark'][r.ds]['A2'][v], 2) for v in 'V0 V1 V2 V3'.split()},
          round(r.F_A1u, 2), {v: round(out['benchmark'][r.ds]['A1u'][v], 2) for v in 'V0 V1 V2 V3'.split()}, flush=True)
# data sets of the study (code of Revision 1)
AF7 = sys.argv[1]; sys.path.insert(0, os.path.join(AF7, 'code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur, gaucher_pipe as gp
def studie(name, X, G, reps, n):
    d = n - 1; o = dict(n=n, p=X.shape[1], G=len(G), A2=vorhersagen(paare(X), d, X.shape[1]), A1u=vorhersagen(paare(X, np.asarray(reps)), d, len(G)))
    out['studien'][name] = o; print(name, n, {v: round(o['A2'][v], 2) for v in 'V0 V1 V2 V3'.split()}, {v: round(o['A1u'][v], 2) for v in 'V0 V1 V2 V3'.split()}, flush=True)
Dm = ng.build('Cleaned_Steps1-4', ycol='H. I. (H. influenzae)'); Xm = Dm['X']; Gm, _ = ao.build_groups(Xm)
print('monoterpenoids keys', list(Dm.keys()), Xm.shape)
studie('Monoterpenoide', Xm, Gm, ur.centrality_reps(Xm, Gm), Xm.shape[0])
Xp = pd.read_csv(os.path.join(AF7, 'revision/lone/desc_pool.csv')); Xp = Xp[[c_ for c_ in Xp.columns if c_ != 'No']].values.astype(float)
Gp, _ = ao.build_groups(Xp); studie('Octapeptide', Xp, Gp, ur.centrality_reps(Xp, Gp), Xp.shape[0])
X, y, mz, _ = gp.load(os.path.join(AF7, 'data/Gaucherdata_paper39_median_normalised.csv')); Xg, _, _ = gp.clean(X, mz)
Gg = gp.build_groups(Xg); Gg = Gg[0] if isinstance(Gg, tuple) else Gg
studie('Gaucher', Xg, Gg, gp.centrality_reps(Xg, Gg), Xg.shape[0])
# sphere model against empirical Gaucher indicator correlations per |rho| class (3000 permutations, run 05.11.2026)
emp = {(0, .1): -0.0131, (.1, .2): -0.0024, (.2, .3): 0.0190, (.3, .5): 0.0641, (.5, .7): 0.1735, (.7, .9): 0.3694, (.9, 1.01): 0.6320}
rg = np.abs(paare(Xg)); tk = tabelle_kugel(38); out['gaucher_klassen'] = []
for (lo, hi), e in emp.items():
    msk = (rg >= lo) & (rg < hi)
    out['gaucher_klassen'].append(dict(klasse=f'{lo}-{hi}', empirisch=e, V0=float(np.mean(k * rg[msk] ** 2)), V2=float(np.mean(np.interp(rg[msk], RG, tk))), V3=float(np.mean(np.interp(rg[msk], RG, TN)))))
    print(out['gaucher_klassen'][-1])
json.dump(out, open(sys.argv[5], 'w'), indent=1)
