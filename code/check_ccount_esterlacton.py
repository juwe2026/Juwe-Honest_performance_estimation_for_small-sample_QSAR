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
import numpy as np, openpyxl
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

wb = openpyxl.load_workbook(_osp.environ.get('QSAR_LITDESK_XLSX', _WORK + '/20260706_31MT_LitDeskriptoren_Analyse.xlsx'), data_only=True)
rows = list(wb["Deskriptoren"].values); hdr = list(rows[0]); data = [r for r in rows[1:] if r[0]]
RED13 = [r[0] for r in list(wb["VIP_reduziert"].values)[1:] if r[0]]
Yc = {'H.I.':'H. I. (H. influenzae)','S.A.':'S. A. (Staph. Aureus)','S.Pneu':'S. Pneu (S. pneumoniae)',
      'S.Pyo':'S. pyo (S. pyogenes)','P.A.':'P. A. (P. aeruginosa)'}
N = len(data)

def col(name): 
    j = hdr.index(name); return np.array([float(r[j]) for r in data])
def yv(s):
    j = hdr.index(Yc[s]); return np.array([1.0 if str(r[j]).strip().lower()=='active' else 0.0 for r in data])

def loo(X, y, nlv):
    yh = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xt = X[tr]; mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd==0] = 1
        m = PLSRegression(n_components=min(nlv, X.shape[1], N-2), scale=False).fit((Xt-mu)/sd, y[tr])
        yh[i] = m.predict((X[i:i+1]-mu)/sd).ravel()[0]
    return yh
def q2(y, yh): return 1 - np.sum((y-yh)**2)/np.sum((y-y.mean())**2)

SETS = {
 "13 (as published)": RED13,
 "12 without C_Count":     [d for d in RED13 if d != "C_Count"],
 "12 without EsterLacton": [d for d in RED13 if d != "EsterLacton_flag"],
 "11 without both":       [d for d in RED13 if d not in ("C_Count","EsterLacton_flag")],
}
for nlv in (1,3):
    print("\n================ nLV = %d ================" % nlv)
    print("%-22s %s" % ("set", "".join("%18s"%s for s in Yc)))
    for name, cols in SETS.items():
        X = np.column_stack([col(c) for c in cols])
        out = []
        for s in Yc:
            y = yv(s); yh = loo(X, y, nlv)
            out.append("Q2%6.3f AUC%5.3f" % (q2(y,yh), roc_auc_score(y,yh)))
        print("%-22s %s" % (name, "".join("%18s"%o for o in out)))
