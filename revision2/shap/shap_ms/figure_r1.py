"""Figure R1 (response to Reviewer 1): SHAP beeswarm plots of the no-selection PLS-DA model (1 LV, all descriptors),
interventional SHAP values (shap.LinearExplainer), monoterpenoids H. influenzae and peptides E. coli. Descriptive, all compounds."""
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
from lade_ms import monoterpene, peptide
from sklearn.cross_decomposition import PLSRegression
import shap
warnings.filterwarnings('ignore')
OUT = sys.argv[1]
fig, axes = plt.subplots(1, 2, figsize=(15, 8.5))
info = []
for ax, (titel, f, ep, lab) in zip(axes, [('a  Monoterpenoids, H. influenzae', monoterpene, 'HI', 'monoterpenoids'), ('b  Peptides, E. coli', peptide, 'EC', 'peptides')]):
    X, Y = f(); y = Y[ep]
    A = X.values; ok = A.std(0, ddof=1) > 1e-6; A = A[:, ok]; names = np.array(X.columns)[ok]
    Z = (A - A.mean(0)) / A.std(0, ddof=1)
    g = routen.gruppierung(Z); grp = {j: k for k, gr in enumerate(g) for j in gr}
    b = PLSRegression(n_components=1, scale=False).fit(Z, y).coef_.ravel()
    phi = shap.LinearExplainer((b, 0.0), shap.maskers.Independent(Z, max_samples=len(Z))).shap_values(Z)
    imp = np.abs(phi).mean(0); top = np.argsort(-imp)[:20]
    order = []; [order.append(grp[j]) for j in top if grp[j] not in order]
    gl = {gk: i + 1 for i, gk in enumerate(order)}
    def kurz(nm):
        src, _, d = nm.partition('::')
        src = {'Mordred': 'Mordred', 'RDKit': 'RDKit', 'PaDEL_1D2D': 'PaDEL', 'PaDEL_3D': 'PaDEL 3D', 'Chemopy_1D2D': 'Chemopy', 'Chemopy_3D': 'Chemopy 3D', 'BlueDesc_mmff94': 'BlueDesc', 'Pybel': 'Pybel'}.get(src, src) if d else ''
        return f"{d} ({src})" if d else nm
    fn = [f"{kurz(names[j])}, group {gl[grp[j]]}" if grp[j] in gl else kurz(names[j]) for j in range(len(names))]
    share = np.sort(imp)[::-1][:20].sum() / imp.sum()
    plt.sca(ax)
    shap.plots.beeswarm(shap.Explanation(values=phi, data=Z, feature_names=fn), max_display=21, ax=ax, show=False, plot_size=None, s=12)
    ax.set_title(f"{titel}\n{Z.shape[1]:,} descriptors, n = {len(y)}\nTop 20: {len(order)} correlation groups, {share*100:.1f} % of total mean |SHAP|", fontsize=10, loc='left')
    ax.set_yticklabels([t.get_text().replace('other features', 'other descriptors').replace('Sum of 2186', 'Sum of 2,186').replace('Sum of 1254', 'Sum of 1,254') for t in ax.get_yticklabels()])
    ax.set_xlabel('SHAP value (contribution to the PLS-DA prediction)', fontsize=10)
    info.append((lab, Z.shape[1], len(order), round(share * 100, 1), float(np.abs(phi[:, top]).max())))
fig.suptitle('Interventional SHAP values of PLS-DA (one latent variable) fitted to all descriptors without selection', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(OUT, dpi=200)
print(info)
