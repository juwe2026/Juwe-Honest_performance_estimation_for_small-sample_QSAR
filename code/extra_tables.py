"""Weitere Route-B-Groessen der Supplementtabellen S7, S8 und S10."""
import sys, json, collections; sys.path.insert(0,'/home/claude/work/recalc')
import numpy as np, core
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

MAXLV = 4
def fit_pred(Xtr, ytr, Xte, nlv):
    mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1
    m = PLSRegression(n_components=min(nlv, Xtr.shape[1], len(ytr)-2), scale=False).fit((Xtr-mu)/sd, ytr)
    return m.predict((Xte-mu)/sd).ravel()

def q2_fixed(X, y, nlv):
    return core.q2(y, core.loo(X, y, nlv))

def doubly_nested(X, y):
    N = len(y); yh = np.zeros(N); chosen = []
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xt, yt = X[tr], y[tr]
        best, bq = 1, -9e9
        for nlv in range(1, MAXLV+1):
            inner = np.zeros(len(yt))
            for j in range(len(yt)):
                t2 = np.ones(len(yt), bool); t2[j] = False
                inner[j] = fit_pred(Xt[t2], yt[t2], Xt[j:j+1], nlv)[0]
            q = core.q2(yt, inner)
            if q > bq: bq, best = q, nlv
        chosen.append(best)
        yh[i] = fit_pred(Xt, yt, X[i:i+1], best)[0]
    c = collections.Counter(chosen)
    return core.q2(y, yh), c.most_common(1)[0][0], min(chosen), max(chosen)

def lpo_c(X, y, nlv=1):
    N = len(y); pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    conc = tie = tot = 0
    for a in pos:
        for b in neg:
            keep = np.ones(N, bool); keep[a] = keep[b] = False
            p = fit_pred(X[keep], y[keep], X[[a, b]], nlv)
            tot += 1
            if p[0] > p[1]: conc += 1
            elif p[0] == p[1]: tie += 1
    return (conc + 0.5*tie)/tot

out = {}
for tag, cols in (("D13", core.RED13), ("D12", core.RED12)):
    X_all = core.X_of(cols); out[tag] = {}
    for s in core.STRAINS:
        y = core.y_of(s); X = X_all
        r = {}
        for nlv in (1, 2, 3):
            r["q2_%dlv" % nlv] = round(float(q2_fixed(X, y, nlv)), 4)
        dq, mode, lo, hi = doubly_nested(X, y)
        r.update(dn_q2=round(float(dq), 4), dn_mode=mode, dn_lo=lo, dn_hi=hi)
        r["auc_loo"] = round(float(roc_auc_score(y, core.loo(X, y, 1))), 4)
        r["lpo_c"] = round(float(lpo_c(X, y, 1)), 4)
        Xs = core.autoscale_all(X)
        m = PLSRegression(n_components=1, scale=False).fit(Xs, y)
        r["r2"] = round(float(1-np.sum((y-m.predict(Xs).ravel())**2)/np.sum((y-y.mean())**2)), 4)
        out[tag][s] = r
        print("%s %-7s 1LV %.2f 2LV %.2f 3LV %.2f | dn %.2f (%d, %d-%d) | AUC %.2f c %.2f R2 %.2f"
              % (tag, s, r["q2_1lv"], r["q2_2lv"], r["q2_3lv"], r["dn_q2"], mode, lo, hi,
                 r["auc_loo"], r["lpo_c"], r["r2"]), flush=True)
json.dump(out, open("recalc/extra_tables.json", "w"), indent=1)
