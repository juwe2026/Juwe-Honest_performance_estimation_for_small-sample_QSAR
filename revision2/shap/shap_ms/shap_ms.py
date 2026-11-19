"""SHAP check on the data sets of the manuscript: PLS-DA, 1 LV, all descriptors (model without selection),
autoscaled (ddof=1 as in common.py, columns with SD < 1e-6 removed), all compounds.
Interventional SHAP values with shap.LinearExplainer (Independent masker), comparison with b_j (x_j - mean_j),
rank correlation with |point-biserial correlation|, redundancy of the top 20 (correlation groups as in the manuscript:
|Pearson| > 0.7 or |Spearman| > 0.7, most-connected-first, routen.gruppierung)."""
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
import sys, json, warnings
import numpy as np
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
import shap
sys.path.insert(0, (_ROOT + '/benchmark')); import routen
from lade_ms import monoterpene, peptide
warnings.filterwarnings('ignore')
res = {}
for name, f in [('Monoterpenoide', monoterpene), ('Peptide', peptide)]:
    X, Y = f(); A = X.values
    sd = A.std(0, ddof=1); ok = sd > 1e-6; A = A[:, ok]
    Z = (A - A.mean(0)) / A.std(0, ddof=1)
    g = routen.gruppierung(Z); grp_of = {j: k for k, grp in enumerate(g) for j in grp}
    C = np.abs(np.nan_to_num(np.corrcoef(Z, rowvar=False)))
    res[name] = dict(n=Z.shape[0], p=Z.shape[1], G=len(g), endpunkte={})
    for ep, y in Y.items():
        pls = PLSRegression(n_components=1, scale=False).fit(Z, y); b = pls.coef_.ravel()
        xty = Z.T @ (y - y.mean())
        r = np.array([np.corrcoef(Z[:, j], y)[0, 1] for j in range(Z.shape[1])])
        ex = shap.LinearExplainer((b, 0.0), shap.maskers.Independent(Z, max_samples=len(Z)))
        phi = ex.shap_values(Z); imp = np.abs(phi).mean(0)
        top = np.argsort(-imp)[:20]; Ct = C[np.ix_(top, top)]
        paare = np.argwhere(np.triu(C > 0.95, 1))
        ratio = np.array([min(imp[i], imp[j]) / max(imp[i], imp[j]) for i, j in paare])
        res[name]['endpunkte'][ep] = dict(aktiv=int(y.sum()),
            kor_coef_XTy=round(float(np.corrcoef(b, xty)[0, 1]), 6),
            max_abw_shap_formel=float(np.abs(phi - b * (Z - Z.mean(0))).max()),
            spearman_imp_r=round(float(spearmanr(imp, np.abs(r)).statistic), 3),
            top20_gruppen=len({grp_of[j] for j in top}),
            top20_paare_r07=int((np.triu(Ct, 1) > 0.7).sum()),
            median_verh_paare_r095=round(float(np.median(ratio)), 3))
        print(name, ep, res[name]['endpunkte'][ep], flush=True)
    print(name, 'n', res[name]['n'], 'p', res[name]['p'], 'G', res[name]['G'], flush=True)
json.dump(res, open(_WORK + '/shap_ms_interventionell.json', 'w'), indent=1)
