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
import sys, warnings, json, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.covariance import LedoitWolf
import shap
warnings.filterwarnings('ignore')
ds = sys.argv[1]; NS = int(sys.argv[2])
A = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_X.csv').values.astype(float)
y = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_y.csv')['y'].values.astype(float)
sd = A.std(0); sd = np.where(sd < 1e-6, 1, sd); Z = (A - A.mean(0)) / sd
r = np.array([np.corrcoef(Z[:, j], y)[0, 1] if Z[:, j].std() > 0 else 0 for j in range(Z.shape[1])])
p = 100; sel = np.argsort(-np.abs(r))[:p]; X = Z[:, sel]
b = PLSRegression(n_components=1, scale=False).fit(X, y).coef_.ravel()
covs = {'stichprobe': np.cov(X, rowvar=False), 'ledoit_wolf': LedoitWolf().fit(X).covariance_}
imp = {}
for name, cov in covs.items():
    for seed in (1, 2):
        np.random.seed(seed); t = time.time()
        ex = shap.LinearExplainer((b, 0.0), (X.mean(0), cov), feature_perturbation='correlation_dependent', nsamples=NS)
        imp[(name, seed)] = np.abs(ex.shap_values(X)).mean(0)
top = lambda v: set(np.argsort(-v)[:10])
sp = lambda a, c: round(float(spearmanr(imp[a], imp[c]).statistic), 3)
ov = lambda a, c: len(top(imp[a]) & top(imp[c]))
o = dict(nsamples=NS, eig_min_stichprobe=float(np.linalg.eigvalsh(covs['stichprobe']).min()), rang_stichprobe=int(np.linalg.matrix_rank(covs['stichprobe'])),
         wiederholung_stichprobe=(sp(('stichprobe',1),('stichprobe',2)), ov(('stichprobe',1),('stichprobe',2))),
         wiederholung_ledoit=(sp(('ledoit_wolf',1),('ledoit_wolf',2)), ov(('ledoit_wolf',1),('ledoit_wolf',2))),
         stichprobe_vs_ledoit=(sp(('stichprobe',1),('ledoit_wolf',1)), ov(('stichprobe',1),('ledoit_wolf',1))))
print(json.dumps(o))
