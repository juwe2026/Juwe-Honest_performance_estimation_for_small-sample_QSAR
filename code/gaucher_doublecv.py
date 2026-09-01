"""Doubly nested analysis of the Smit 2007 Gaucher data (Table S12, last row).

The main Gaucher results fix the number of latent variables a priori. This script instead
lets an inner cross-validation choose it, so that the component number is never selected on
data used to score the model:

  outer loop   leave-one-out over all 39 samples;
  selection    the route's descriptor selection is run on the 38 training samples of the
               outer fold, exactly as in the fixed-component analysis;
  inner loop   a second leave-one-out over those 38 samples, run once for every candidate
               component number from 1 to MAXLV on the already-selected descriptors;
  choice       the candidate with the highest inner Q2 wins, ties going to the smaller
               number, and the outer prediction is made with it.

The unsupervised correlation grouping is the one fixed on all 39 samples, as everywhere else
in this study; only the label-dependent steps are recomputed inside every fold.

Usage:  python code/gaucher_doublecv.py [ROUTE ...]     (default: A1 A1u A2)
Output: results/recalc/gaucher_doublecv.json and a printed summary.

This reproduces the archived results/raw_permutations/gaucher_doublecv.pkl exactly, including
the per-fold record of which component number was chosen.
"""
import os as _os
import sys, json, time
import numpy as np

_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
sys.path.insert(0, _HERE)

import gaucher_pipe as gp
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

CSV = _os.environ.get(
    "QSAR_GAUCHER_CSV",
    _os.path.join(_DATA, "Gaucherdata_paper39_median_normalised.csv"),
)
MAXLV = 5
ROUTES = sys.argv[1:] or ["A1", "A1u", "A2"]


def fit_pred(Xtr, ytr, Xte, nlv):
    mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
    m = PLSRegression(n_components=min(nlv, Xtr.shape[1], len(ytr) - 2),
                      scale=False).fit((Xtr - mu) / sd, ytr)
    return m.predict((Xte - mu) / sd).ravel()


def doubly_nested(Xc, y, route, groups, reps_u):
    n = len(y); yh = np.zeros(n); chosen = []
    for i in range(n):
        tr = np.ones(n, bool); tr[i] = False
        Xt, yt = Xc[tr], y[tr]
        cols, _ = gp.select(route, Xt, yt, groups, reps_u)
        A = Xt[:, cols]
        best, bq = 1, -9e9
        for nlv in range(1, MAXLV + 1):
            inner = np.zeros(len(yt))
            for j in range(len(yt)):
                t2 = np.ones(len(yt), bool); t2[j] = False
                inner[j] = fit_pred(A[t2], yt[t2], A[j:j + 1], nlv)[0]
            q = gp.q2(yt, inner)
            if q > bq:                      # strict: ties go to the smaller number
                bq, best = q, nlv
        chosen.append(best)
        yh[i] = fit_pred(A, yt, Xc[i:i + 1, cols], best)[0]
    return gp.q2(y, yh), roc_auc_score(y, yh), chosen, yh


def main():
    t0 = time.time()
    X, y, mz, _ = gp.load(CSV)
    Xc, mzc, _ = gp.clean(X, mz)
    groups, _ = gp.build_groups(Xc)
    reps_u = gp.centrality_reps(Xc, groups)
    print("%d samples, %d m/z, %d cases, %d correlation groups, candidates 1-%d [%.0f s]"
          % (Xc.shape[0], Xc.shape[1], int(y.sum()), len(groups), MAXLV, time.time() - t0),
          flush=True)

    out = {"maxlv": MAXLV, "routes": {}}
    for route in ROUTES:
        t = time.time()
        q, a, chosen, yh = doubly_nested(Xc, y, route, groups, reps_u)
        m = gp.metrics(y, yh)
        out["routes"][route] = {"q2": float(q), "auc": float(a),
                                "chosen": [int(c) for c in chosen], "metrics": m}
        vals, cnt = np.unique(chosen, return_counts=True)
        print("%-4s doubly nested Q2 %+.4f  AUC %.4f  components chosen: %s  [%.0f s]"
              % (route, q, a, dict(zip(vals.tolist(), cnt.tolist())), time.time() - t),
              flush=True)

    path = _os.path.join(_RES, "gaucher_doublecv.json")
    json.dump(out, open(path, "w"), indent=1, default=float)
    print("done after %.0f s -> %s" % (time.time() - t0, path))


if __name__ == "__main__":
    main()
