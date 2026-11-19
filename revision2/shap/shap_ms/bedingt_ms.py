"""Conditional SHAP values on the data of the manuscript (subsets of the descriptors): dependence on the covariance model
and attribution to unused descriptors. Same procedure as ../shap_pruef/stabil.py and nur5.py."""
import sys, json, warnings
import numpy as np
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.covariance import LedoitWolf
import shap
from lade_ms import monoterpene, peptide, _WORK
warnings.filterwarnings('ignore')
name, ep, NS = sys.argv[1], sys.argv[2], int(sys.argv[3])
X, Y = {'Monoterpenoide': monoterpene, 'Peptide': peptide}[name](); y = Y[ep]
A = X.values; ok = A.std(0, ddof=1) > 1e-6; A = A[:, ok]; Z = (A - A.mean(0)) / A.std(0, ddof=1)
r = np.array([np.corrcoef(Z[:, j], y)[0, 1] for j in range(Z.shape[1])])
o = dict(datensatz=name, endpunkt=ep, n=len(y), nsamples=NS)
# a) p = 100, two covariance models, two seeds each
sel = np.argsort(-np.abs(r))[:100]; Xs = Z[:, sel]
b = PLSRegression(n_components=1, scale=False).fit(Xs, y).coef_.ravel()
covs = {'stichprobe': np.cov(Xs, rowvar=False), 'ledoit_wolf': LedoitWolf().fit(Xs).covariance_}
o['rang_stichprobenkovarianz'] = int(np.linalg.matrix_rank(covs['stichprobe']))
imp = {}
for cn, cov in covs.items():
    for seed in (1, 2):
        np.random.seed(seed)
        ex = shap.LinearExplainer((b, 0.0), (Xs.mean(0), cov), feature_perturbation='correlation_dependent', nsamples=NS)
        imp[(cn, seed)] = np.abs(ex.shap_values(Xs)).mean(0)
top = lambda v: set(np.argsort(-v)[:10])
cmp_ = lambda a, c: (round(float(spearmanr(imp[a], imp[c]).statistic), 3), len(top(imp[a]) & top(imp[c])))
o['wiederholung_stichprobe'] = cmp_(('stichprobe', 1), ('stichprobe', 2))
o['wiederholung_ledoit'] = cmp_(('ledoit_wolf', 1), ('ledoit_wolf', 2))
o['stichprobe_vs_ledoit'] = cmp_(('stichprobe', 1), ('ledoit_wolf', 1))
# b) model uses 5 of 50 descriptors
sel = np.argsort(-np.abs(r))[:50]; Xs = Z[:, sel]; use = np.arange(0, 50, 10); unused = np.setdiff1d(np.arange(50), use)
bb = np.zeros(50); bb[use] = PLSRegression(n_components=1, scale=False).fit(Xs[:, use], y).coef_.ravel()
for cn, cov in [('ledoit_wolf', LedoitWolf().fit(Xs).covariance_), ('stichprobe', np.cov(Xs, rowvar=False))]:
    for seed in (1, 2):
        np.random.seed(seed)
        ex = shap.LinearExplainer((bb, 0.0), (Xs.mean(0), cov), feature_perturbation='correlation_dependent', nsamples=10000)
        phi = np.abs(ex.shap_values(Xs)).mean(0)
        o[f'nur5_{cn}_seed{seed}'] = (round(float(phi[unused].sum() / phi.sum()), 3), int((phi[unused] > phi[use].min()).sum()))
print(json.dumps(o)); json.dump(o, open(f'{_WORK}/bedingt_{name}_{ep}.json', 'w'), indent=1)
