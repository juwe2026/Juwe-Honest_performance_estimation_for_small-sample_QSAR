"""VIP-Werte fuer Fig. S2b und Fig. S3a: apparent und Mittel ueber die LOO-Folds,
zwei latente Variablen, wie in 20260710_VIP_EffectSize_Tabellen.xlsx."""
import sys, json; sys.path.insert(0, '/home/claude/work/recalc')
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
    wb = openpyxl.load_workbook("/mnt/user-data/uploads/QSAR_Monoterpene_Gesamtsicherung_2026-08-11/20260710_VIP_EffectSize_Tabellen.xlsx", data_only=True)
    v13 = vip_set(core.RED13)
    err = 0.0
    for s in core.STRAINS:
        for r in list(wb["VIP_Lit_"+s].values)[1:]:
            if r[0] in core.RED13 and r[1] is not None:
                err += abs(v13[s][r[0]]["apparent"]-r[1]) + abs(v13[s][r[0]]["nested"]-r[2])
    print("Kontrolle 13 Deskriptoren, Summe |Abweichung| zur publizierten Tabelle: %.4f" % err)
    v12 = vip_set(core.RED12)
    json.dump({"nlv": NLV, "D13": v13, "D12": v12},
              open("figs/data_recalc/vip_routeB.json", "w"), indent=1)
    print("geschrieben: figs/data_recalc/vip_routeB.json")
