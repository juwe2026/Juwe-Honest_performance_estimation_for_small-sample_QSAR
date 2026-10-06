# Revision (EW 03.10.2026): Zahl der latenten Variablen fuer die Arme ohne Selektion (Table S5).
# Apparent (Anpassung und Vorhersage an allen 31 Verbindungen), genestet (LOO; ohne Selektion bereits ehrlich)
# fuer 1 bis 4 LV, und doppelt genestet: Komponentenzahl 1-4 per innerer LOO-CV auf den Trainingsverbindungen
# jedes aeusseren Folds gewaehlt. Keine Permutationen, deterministisch.
exec(open('/home/claude/rev/common.py').read())
from sklearn.metrics import roc_auc_score
import json, collections
def fit_pred(Xtr, ytr, Xte, L):
    mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); ok = sd > 1e-6
    m = PLSRegression(min(L, len(ytr) - 2), scale=False).fit((Xtr[:, ok] - mu[ok]) / sd[ok], ytr)
    return m.predict((Xte[:, ok] - mu[ok]) / sd[ok]).ravel()
out = {}
for nm, Xm in (('all2208', X), ('reps182', X[:, REPS_U])):
    out[nm] = {}
    for k in KEYS:
        y = Y[k]; r = {'apparent': {}, 'nested': {}}
        for L in (1, 2, 3, 4):
            ya = fit_pred(Xm, y, Xm, L); yh = loo_pls(Xm, y, L)
            r['apparent'][L] = dict(q2=q2(y, ya), auc=roc_auc_score(y, ya))
            r['nested'][L] = dict(q2=q2(y, yh), auc=roc_auc_score(y, yh))
        yh = np.zeros(N); ch = []
        for i in range(N):
            tr = np.ones(N, bool); tr[i] = False; Xt, yt = Xm[tr], y[tr]; best = None
            for L in (1, 2, 3, 4):
                p = np.array([fit_pred(np.delete(Xt, j, 0), np.delete(yt, j), Xt[j:j + 1], L)[0] for j in range(len(yt))])
                q = q2(yt, p)
                if best is None or q > best[1]: best = (L, q)
            ch.append(best[0]); yh[i] = fit_pred(Xt, yt, Xm[i:i + 1], best[0])[0]
        c = collections.Counter(ch)
        r['doubly'] = dict(q2=q2(y, yh), auc=roc_auc_score(y, yh), choice={str(a): b for a, b in sorted(c.items())},
                           median=int(np.median(ch)), lo=int(min(ch)), hi=int(max(ch)))
        out[nm][k] = r
        print(nm, k, 'app', [round(r['apparent'][L]['q2'], 2) for L in (1, 2, 3, 4)], 'nest', [round(r['nested'][L]['q2'], 2) for L in (1, 2, 3, 4)],
              'doppelt %.3f' % r['doubly']['q2'], dict(c), flush=True)
json.dump(out, open('/home/claude/rev/out/lv_noselection.json', 'w'), indent=1, default=float)
