"""Figure R1, version with two x scales: at the top the 20 descriptors with the largest mean |SHAP| (own, fine scale),
at the bottom the sum of all other descriptors (own scale). Otherwise as figure_r1.py."""
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
from matplotlib.gridspec import GridSpec
sys.path.insert(0, (_ROOT + '/benchmark')); import routen
sys.path.insert(0, _ROOT + '/revision2/shap/shap_ms')
from lade_ms import monoterpene, peptide
from sklearn.cross_decomposition import PLSRegression
import shap
warnings.filterwarnings('ignore')
OUT = sys.argv[1]
SRC = {'Mordred': 'Mordred', 'RDKit': 'RDKit', 'PaDEL_1D2D': 'PaDEL', 'PaDEL_3D': 'PaDEL 3D', 'Chemopy_1D2D': 'Chemopy', 'Chemopy_3D': 'Chemopy 3D', 'BlueDesc_mmff94': 'BlueDesc', 'Pybel': 'Pybel'}
def kurz(nm):
    s, _, d = nm.partition('::'); return f"{d} ({SRC.get(s, s)})" if d else nm
fig = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 2, height_ratios=[20, 2.4], hspace=0.32, wspace=0.75, figure=fig)
from matplotlib.ticker import MaxNLocator
AXES = []
info = []
for col, (titel, f, ep) in enumerate([('a  Monoterpenoids, H. influenzae', monoterpene, 'HI'), ('b  Peptides, E. coli', peptide, 'EC')]):
    X, Y = f(); y = Y[ep]
    A = X.values; ok = A.std(0, ddof=1) > 1e-6; A = A[:, ok]; names = np.array(X.columns)[ok]
    Z = (A - A.mean(0)) / A.std(0, ddof=1)
    g = routen.gruppierung(Z); grp = {j: k for k, gr in enumerate(g) for j in gr}
    b = PLSRegression(n_components=1, scale=False).fit(Z, y).coef_.ravel()
    phi = shap.LinearExplainer((b, 0.0), shap.maskers.Independent(Z, max_samples=len(Z))).shap_values(Z)
    imp = np.abs(phi).mean(0); top = np.argsort(-imp)[:20]; rest = np.setdiff1d(np.arange(Z.shape[1]), top)
    order = []; [order.append(grp[j]) for j in top if grp[j] not in order]
    gl = {gk: i + 1 for i, gk in enumerate(order)}
    share = imp[top].sum() / imp.sum()
    ax1 = fig.add_subplot(gs[0, col]); ax2 = fig.add_subplot(gs[1, col])
    e1 = shap.Explanation(values=phi[:, top], data=Z[:, top], feature_names=[f"{kurz(names[j])}, group {gl[grp[j]]}" for j in top])
    shap.plots.beeswarm(e1, max_display=20, ax=ax1, show=False, plot_size=None, s=12, group_remaining_features=False)
    lim1 = np.abs(phi[:, top]).max() * 1.15; ax1.set_xlim(-lim1, lim1)
    ax1.xaxis.set_major_locator(MaxNLocator(nbins=5, symmetric=True))
    ax1.set_xlabel('SHAP value of the single descriptor (fine scale)', fontsize=9)
    ax1.set_title(f"{titel}\n{Z.shape[1]:,} descriptors, n = {len(y)}\nTop 20: {len(order)} correlation groups, {share*100:.1f} % of total mean |SHAP|", fontsize=10, loc='left')
    e2 = shap.Explanation(values=phi[:, rest].sum(1, keepdims=True), data=None, feature_names=[f"Sum of the other {len(rest):,} descriptors"])
    shap.plots.beeswarm(e2, max_display=1, ax=ax2, show=False, plot_size=None, s=12, color_bar=False, color='#7f7f7f')
    lim2 = np.abs(phi[:, rest].sum(1)).max() * 1.1; ax2.set_xlim(-lim2, lim2)
    ax2.xaxis.set_major_locator(MaxNLocator(nbins=7, symmetric=True))
    ax2.set_xlabel(f'SHAP value, sum of the other descriptors (scale about {round(lim2 / lim1, -1):.0f} times coarser)', fontsize=9)
    AXES.append((ax1, ax2))
    info.append((titel, round(lim1, 4), round(lim2, 3), round(lim2 / lim1)))
fig.canvas.draw()
for a1, a2 in AXES:
    p1 = a1.get_position(); p2 = a2.get_position(); a2.set_position([p1.x0, p2.y0, p1.width, p2.height])
fig.suptitle('Interventional SHAP values of PLS-DA (one latent variable) fitted to all descriptors without selection', fontsize=12, y=0.995)
fig.savefig(OUT, dpi=200, bbox_inches='tight'); print(info)
