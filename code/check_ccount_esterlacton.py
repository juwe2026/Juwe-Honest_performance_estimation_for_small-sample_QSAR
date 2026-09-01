import numpy as np, openpyxl
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

wb = openpyxl.load_workbook("/mnt/user-data/uploads/QSAR_Monoterpene_Gesamtsicherung_2026-08-11/20260706_31MT_LitDeskriptoren_Analyse.xlsx", data_only=True)
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
 "13 (wie publiziert)": RED13,
 "12 ohne C_Count":     [d for d in RED13 if d != "C_Count"],
 "12 ohne EsterLacton": [d for d in RED13 if d != "EsterLacton_flag"],
 "11 ohne beide":       [d for d in RED13 if d not in ("C_Count","EsterLacton_flag")],
}
for nlv in (1,3):
    print("\n================ nLV = %d ================" % nlv)
    print("%-22s %s" % ("Satz", "".join("%18s"%s for s in Yc)))
    for name, cols in SETS.items():
        X = np.column_stack([col(c) for c in cols])
        out = []
        for s in Yc:
            y = yv(s); yh = loo(X, y, nlv)
            out.append("Q2%6.3f AUC%5.3f" % (q2(y,yh), roc_auc_score(y,yh)))
        print("%-22s %s" % (name, "".join("%18s"%o for o in out)))
