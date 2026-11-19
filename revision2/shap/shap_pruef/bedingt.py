"""Conditional (correlation_dependent) SHAP values for PLS-DA with 1 LV: computational cost and dependence on the
model of the dependence (covariance estimate). Subsets of the descriptors, because the full set cannot be computed."""
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
import json, sys, time, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.covariance import LedoitWolf
import shap
warnings.filterwarnings('ignore')
ds = sys.argv[1]
A = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_X.csv').values.astype(float)
y = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_y.csv')['y'].values.astype(float)
mu = A.mean(0); sd = A.std(0); sd = np.where(sd < 1e-6, 1, sd); Z = (A - mu) / sd
r = np.array([np.corrcoef(Z[:, j], y)[0, 1] if Z[:, j].std() > 0 else 0 for j in range(Z.shape[1])])
out = {}
# 1) computing time as a function of p
for p in (25, 50, 100):
    sel = np.argsort(-np.abs(r))[:p]
    pls = PLSRegression(n_components=1, scale=False).fit(Z[:, sel], y); b = pls.coef_.ravel()
    cov = np.cov(Z[:, sel], rowvar=False)
    t = time.time()
    ex = shap.LinearExplainer((b, 0.0), (np.zeros(p), cov), feature_perturbation='correlation_dependent', nsamples=50)
    out[f'zeit_p{p}_50samples_s'] = round(time.time() - t, 2)
    print(p, out[f'zeit_p{p}_50samples_s'], flush=True)
# 2) same data, two models of the dependence (sample covariance + 1e-6 as in shap, Ledoit-Wolf)
p = 100; sel = np.argsort(-np.abs(r))[:p]
pls = PLSRegression(n_components=1, scale=False).fit(Z[:, sel], y); b = pls.coef_.ravel()
imp = {}
for name, cov in [('stichprobe', np.cov(Z[:, sel], rowvar=False)), ('ledoit_wolf', LedoitWolf().fit(Z[:, sel]).covariance_)]:
    ex = shap.LinearExplainer((b, 0.0), (Z[:, sel].mean(0), cov), feature_perturbation='correlation_dependent', nsamples=1000)
    phi = ex.shap_values(Z[:, sel]); imp[name] = np.abs(phi).mean(0)
imp['interventionell'] = np.abs(b * (Z[:, sel] - Z[:, sel].mean(0))).mean(0)
k = 10; tops = {n: set(np.argsort(-v)[:k]) for n, v in imp.items()}
out['p100_spearman_stichprobe_ledoit'] = float(spearmanr(imp['stichprobe'], imp['ledoit_wolf']).statistic)
out['p100_spearman_stichprobe_interv'] = float(spearmanr(imp['stichprobe'], imp['interventionell']).statistic)
out['p100_spearman_ledoit_interv'] = float(spearmanr(imp['ledoit_wolf'], imp['interventionell']).statistic)
out['p100_top10_ueberlapp_stichprobe_ledoit'] = len(tops['stichprobe'] & tops['ledoit_wolf'])
out['p100_top10_ueberlapp_stichprobe_interv'] = len(tops['stichprobe'] & tops['interventionell'])
out['p100_top10_ueberlapp_ledoit_interv'] = len(tops['ledoit_wolf'] & tops['interventionell'])
# 3) model uses only 5 descriptors, conditional SHAP values over 50 descriptors (45 not used by the model)
p = 50; sel = np.argsort(-np.abs(r))[:p]
use = np.arange(0, 50, 10)          # 5 descriptors used, spread over the ranking
pls = PLSRegression(n_components=1, scale=False).fit(Z[:, sel[use]], y)
b = np.zeros(p); b[use] = pls.coef_.ravel()
cov = LedoitWolf().fit(Z[:, sel]).covariance_
ex = shap.LinearExplainer((b, 0.0), (Z[:, sel].mean(0), cov), feature_perturbation='correlation_dependent', nsamples=1000)
phi = np.abs(ex.shap_values(Z[:, sel])).mean(0)
unused = np.setdiff1d(np.arange(p), use)
out['nur5_genutzt_anteil_wichtigkeit_auf_ungenutzten'] = float(phi[unused].sum() / phi.sum())
out['nur5_genutzt_ungenutzte_mit_wichtigkeit_ueber_kleinster_genutzter'] = int((phi[unused] > phi[use].min()).sum())
phi_i = np.abs(b * (Z[:, sel] - Z[:, sel].mean(0))).mean(0)
out['nur5_interventionell_wichtigkeit_ungenutzte'] = float(phi_i[unused].sum())
print(json.dumps(out, indent=1))
json.dump(out, open(f'{_WORK}/bedingt_{ds}.json', 'w'), indent=1)
