"""VIP values for Fig. S2b and Fig. S3a: apparent and mean over the LOO folds,
two latent variables, as in 20260710_VIP_EffectSize_Tabellen.xlsx."""
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
import sys, json; sys.path.insert(0, (_ROOT + '/code'))
import numpy as np, openpyxl, core
from sklearn.cross_decomposition import PLSRegression

NLV = 2
def vip_set(cols):
    X = core.X_of(cols); out = {}
    for s in core.STRAINS:
        y = core.y_of(s)
        Xs = core.autoscale_all(X)
        app = core.vip(PLSRegression(n_components=NLV, scale=False).fit(Xs, y), Xs)
        acc = np.zeros(len(cols))
        for i in range(core.N):
            tr = np.ones(core.N, bool); tr[i] = False
            Xt = X[tr]; mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd == 0] = 1
            Z = (Xt-mu)/sd
            acc += core.vip(PLSRegression(n_components=NLV, scale=False).fit(Z, y[tr]), Z)
        out[s] = {c: dict(apparent=float(app[k]), nested=float(acc[k]/core.N))
                  for k, c in enumerate(cols)}
    return out

if __name__ == "__main__":
    wb = openpyxl.load_workbook(_osp.environ.get('QSAR_VIP_XLSX', _ROOT + '/data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx'), data_only=True)
    v13 = vip_set(core.RED13)
    err = 0.0
    for s in core.STRAINS:
        for r in list(wb["VIP_Lit_"+s].values)[1:]:
            if r[0] in core.RED13 and r[1] is not None:
                err += abs(v13[s][r[0]]["apparent"]-r[1]) + abs(v13[s][r[0]]["nested"]-r[2])
    print("Check 13 descriptors, sum |deviation| from the published table: %.4f" % err)
    v12 = vip_set(core.RED12)
    json.dump({"nlv": NLV, "D13": v13, "D12": v12},
              open(_WORK + '/vip_routeB.json', 'w'), indent=1)
    print("written: figs/data_recalc/vip_routeB.json")
