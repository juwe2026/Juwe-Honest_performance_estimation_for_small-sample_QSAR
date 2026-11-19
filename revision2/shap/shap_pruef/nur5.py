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
import sys, warnings, json
import numpy as np, pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.covariance import LedoitWolf
import shap
warnings.filterwarnings('ignore')
ds = sys.argv[1]; NS = int(sys.argv[2])
A = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_X.csv').values.astype(float)
y = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_y.csv')['y'].values.astype(float)
sd = A.std(0); sd = np.where(sd < 1e-6, 1, sd); Z = (A - A.mean(0)) / sd
r = np.array([np.corrcoef(Z[:, j], y)[0, 1] if Z[:, j].std() > 0 else 0 for j in range(Z.shape[1])])
p = 50; sel = np.argsort(-np.abs(r))[:p]; X = Z[:, sel]
use = np.arange(0, 50, 10); unused = np.setdiff1d(np.arange(p), use)
b = np.zeros(p); b[use] = PLSRegression(n_components=1, scale=False).fit(X[:, use], y).coef_.ravel()
o = {}
for name, cov in [('ledoit_wolf', LedoitWolf().fit(X).covariance_), ('stichprobe', np.cov(X, rowvar=False))]:
    for seed in (1, 2):
        np.random.seed(seed)
        ex = shap.LinearExplainer((b, 0.0), (X.mean(0), cov), feature_perturbation='correlation_dependent', nsamples=NS)
        phi = np.abs(ex.shap_values(X)).mean(0)
        o[f'{name}_seed{seed}'] = dict(anteil_ungenutzt=round(float(phi[unused].sum() / phi.sum()), 3),
                                       ungenutzte_ueber_kleinstem_genutzten=int((phi[unused] > phi[use].min()).sum()))
print(json.dumps(dict(nsamples=NS, **o)))
