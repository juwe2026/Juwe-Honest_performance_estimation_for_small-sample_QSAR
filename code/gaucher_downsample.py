"""Detection-limit analysis on the Smit 2007 Gaucher data (Table S13).

From the 39 samples (20 controls, 19 Gaucher cases), k of the 19 cases are drawn at
random without replacement while all 20 controls are retained, giving a subsample of
20 + k samples whose active fraction matches one of the monoterpenoid endpoints. Route
A1 with one latent variable is then run on the subsample under fully nested
leave-one-out validation, exactly as in the main analysis. The unsupervised correlation
grouping is the one fixed on all 39 samples, as everywhere else in this study; only the
label-dependent steps are recomputed inside every fold.

Usage:  python code/gaucher_downsample.py [NDRAWS]        (default 2000)
Output: results/recalc/gaucher_downsample.json  and a printed summary table.

Reproducibility: draw d for case count k uses numpy.random.default_rng([SEED, k, d]),
so every subsample is determined by (SEED, k, d) alone and any single draw can be
regenerated in isolation.
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
from sklearn.metrics import roc_auc_score

CSV = _os.environ.get(
    "QSAR_GAUCHER_CSV",
    _os.path.join(_DATA, "Gaucherdata_paper39_median_normalised.csv"),
)
SEED = 20260829
KS = [11, 10, 6, 5, 4]          # matched to S. pyogenes, H. influenzae, S. pneumoniae, S. aureus
NDRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
ROUTE, NLV = "A1", 1


def main():
    t0 = time.time()
    X, y, mz, _ = gp.load(CSV)
    Xc, mzc, _ = gp.clean(X, mz)
    groups, _ = gp.build_groups(Xc)
    reps_u = gp.centrality_reps(Xc, groups)
    cases = np.where(y == 1)[0]
    ctrl = np.where(y == 0)[0]
    print("%d samples, %d m/z, %d cases, %d controls, %d correlation groups [%.0f s]"
          % (Xc.shape[0], Xc.shape[1], len(cases), len(ctrl), len(groups), time.time() - t0),
          flush=True)

    def nested(Xs, ys):
        yh = gp.loo_nested(Xs, ys, ROUTE, groups, reps_u, nlv=NLV)
        return gp.q2(ys, yh), roc_auc_score(ys, yh)

    out = {"seed": SEED, "ndraws": NDRAWS, "route": ROUTE, "nlv": NLV, "k": {}}

    q, a = nested(Xc, y)                      # k = 19: the full data, a single value
    out["k"]["19"] = {"n_draws": 1, "q2": [float(q)], "auc": [float(a)]}
    print("k=19  n=39  Q2 %+.4f  AUC %.4f  (single value, no sampling)" % (q, a), flush=True)

    for k in KS:
        t = time.time()
        q2s = np.empty(NDRAWS)
        aucs = np.empty(NDRAWS)
        for d in range(NDRAWS):
            rng = np.random.default_rng([SEED, k, d])
            sel = np.sort(np.concatenate([rng.choice(cases, size=k, replace=False), ctrl]))
            q2s[d], aucs[d] = nested(Xc[sel], y[sel])
        out["k"][str(k)] = {"n_draws": NDRAWS, "q2": q2s.tolist(), "auc": aucs.tolist()}
        print("k=%-2d  n=%d  Q2 median %+.4f  (P5 %+.4f, P95 %+.4f)  AUC median %.4f  "
              "Q2>0 %5.1f%%  Q2>0.30 %5.1f%%  [%.0f s]"
              % (k, 20 + k, np.median(q2s), np.percentile(q2s, 5), np.percentile(q2s, 95),
                 np.median(aucs), 100 * np.mean(q2s > 0), 100 * np.mean(q2s > 0.30),
                 time.time() - t), flush=True)

    path = _os.path.join(_RES, "gaucher_downsample.json")
    json.dump(out, open(path, "w"))
    print("done after %.0f s -> %s" % (time.time() - t0, path))


if __name__ == "__main__":
    main()
