"""Beeswarm plot (interventional SHAP values) for PLS-DA, 1 LV, all descriptors, monoterpenoids H. influenzae.
For viewing only, descriptive, model on all 31 compounds."""
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
import sys, warnings, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, (_ROOT + '/benchmark')); import routen
from lade_ms import monoterpene
from sklearn.cross_decomposition import PLSRegression
import shap
warnings.filterwarnings('ignore')
ep = sys.argv[1]; out = sys.argv[2]
X, Y = monoterpene(); y = Y[ep]
A = X.values; ok = A.std(0, ddof=1) > 1e-6; A = A[:, ok]; names = np.array(X.columns)[ok]
Z = (A - A.mean(0)) / A.std(0, ddof=1)
g = routen.gruppierung(Z); grp = {j: k for k, gr in enumerate(g) for j in gr}
b = PLSRegression(n_components=1, scale=False).fit(Z, y).coef_.ravel()
phi = shap.LinearExplainer((b, 0.0), shap.maskers.Independent(Z, max_samples=len(Z))).shap_values(Z)
imp = np.abs(phi).mean(0); top = np.argsort(-imp)[:20]
# number the groups of the top 20 in the order of their first occurrence
order = []; [order.append(grp[j]) for j in top if grp[j] not in order]
lab = {gk: i + 1 for i, gk in enumerate(order)}
fn = [f"{names[j].split('::')[-1]}  [group {lab[grp[j]]}]" if grp[j] in lab else names[j].split('::')[-1] for j in range(len(names))]
ex = shap.Explanation(values=phi, data=Z, feature_names=fn)
plt.figure()
shap.plots.beeswarm(ex, max_display=21, show=False)
plt.title(f'Monoterpenoids, {ep}: PLS-DA (1 LV, 2,206 descriptors), interventional SHAP values\nthe top 20 come from {len(order)} correlation groups', fontsize=10)
plt.gcf().set_size_inches(10, 9); plt.tight_layout(); plt.savefig(out, dpi=150)
print('groups in the top 20:', len(order), [lab[grp[j]] for j in top])
