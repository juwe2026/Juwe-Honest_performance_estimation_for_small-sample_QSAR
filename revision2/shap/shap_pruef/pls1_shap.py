"""Check: what does SHAP give for PLS-DA with one latent variable on all descriptors (model 'ohneAuswahl', no selection)?
Data: ChEMBL benchmark, datasets_k (main analysis). No result enters the hypothesis tests.
"""
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
from scipy.stats import spearmanr, pearsonr
from sklearn.cross_decomposition import PLSRegression
import shap
warnings.filterwarnings('ignore')

def lade(ds):
    X = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_X.csv'); cols = list(X.columns)
    y = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_y.csv')['y'].values.astype(float)
    info = json.load(open(f'{_ROOT}/benchmark/datasets_k/{ds}_info.json'))
    g = [[cols.index(c) for c in grp if c in cols] for grp in info['gruppen']]
    return X.values.astype(float), y, g

res = {}
for ds in sys.argv[1:]:
    A, y, g = lade(ds)
    mu = A.mean(0); sd = A.std(0); sd = np.where(sd < 1e-6, 1, sd); Z = (A - mu) / sd
    pls = PLSRegression(n_components=1, scale=False).fit(Z, y)
    b = pls.coef_.ravel()
    xty = Z.T @ (y - y.mean())
    r_pb = np.array([pearsonr(Z[:, j], y)[0] if Z[:, j].std() > 0 else 0 for j in range(Z.shape[1])])
    # interventional SHAP values (shap.LinearExplainer, Independent masker)
    ex = shap.LinearExplainer((b, float(pls.intercept_.ravel()[0]) if np.ndim(pls.intercept_) else float(pls.intercept_)), shap.maskers.Independent(Z, max_samples=len(Z)))
    phi = ex.shap_values(Z)
    phi_formel = b * (Z - Z.mean(0))
    imp = np.abs(phi).mean(0)
    # groups with at least 2 descriptors: spread of the importance within the group
    C = np.corrcoef(Z, rowvar=False)
    paare = np.argwhere(np.triu(np.abs(np.nan_to_num(C)) > 0.95, 1))
    ratio = np.array([min(imp[i], imp[j]) / max(imp[i], imp[j]) for i, j in paare if max(imp[i], imp[j]) > 0])
    top = np.argsort(-imp)[:20]
    Ct = np.abs(np.nan_to_num(C[np.ix_(top, top)])); np.fill_diagonal(Ct, 0)
    # groups of the top 20 (correlation groups from the specification)
    grp_of = {j: k for k, grp in enumerate(g) for j in grp}
    res[ds] = dict(n=len(y), p=Z.shape[1], G=len(g),
                   kor_coef_XTy=float(np.corrcoef(b, xty)[0, 1]),
                   max_abw_shap_formel=float(np.abs(phi - phi_formel).max()),
                   spearman_imp_rpb=float(spearmanr(imp, np.abs(r_pb)).statistic),
                   paare_r095=int(len(paare)), median_verhaeltnis_paar=float(np.median(ratio)) if len(ratio) else None,
                   anteil_paare_verh_ueber_0_8=float((ratio > 0.8).mean()) if len(ratio) else None,
                   top20_gruppen=len({grp_of.get(j) for j in top}),
                   top20_paare_r07=int((np.triu(Ct, 1) > 0.7).sum()))
    print(ds, res[ds], flush=True)
json.dump(res, open(_WORK + '/pls1_shap_ergebnis.json', 'w'), indent=1)
