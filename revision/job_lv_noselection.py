# Revision (EW 03.10.2026): number of latent variables for the arms without selection (Table S5 of Revision 1).
# Apparent (fit and prediction on all 31 compounds), nested (LOO, already honest without selection)
# for 1 to 4 LV, and double nested: number of components 1-4 chosen by inner LOO-CV on the training compounds
# of every outer fold. No permutations, deterministic.
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
exec(open((_REV + '/common.py')).read())
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
              'doubly nested %.3f' % r['doubly']['q2'], dict(c), flush=True)
json.dump(out, open((_OUT + '/lv_noselection.json'), 'w'), indent=1, default=float)
