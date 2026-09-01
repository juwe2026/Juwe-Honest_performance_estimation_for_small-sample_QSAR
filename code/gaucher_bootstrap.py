"""Why the Gaucher subsamples are drawn without replacement (Table S13, legend).

Only 20 controls exist, so a subsample cannot match both the number of actives and the
sample size of our monoterpenoid endpoints at once. Drawing controls with replacement would
remove that constraint, and this script measures what it would cost.

For each case count k the controls are bootstrapped to bring the subsample to n = 31: k of
the 19 patients are drawn without replacement, and 31 - k controls are drawn from the 20 with
replacement. Route A1 with one latent variable is then run under fully nested leave-one-out
validation, exactly as in the without-replacement analysis of `gaucher_downsample.py`. Because
a duplicated control's identical twin stays in the training set whenever leave-one-out holds
it out, the model has already seen the held-out sample, and the estimate is inflated. The
script also records how many of the drawn controls are duplicates.

Usage:  python code/gaucher_bootstrap.py [NDRAWS]       (default 2000)
Output: results/recalc/gaucher_bootstrap.json and a printed summary.

Reproducibility: draw d for case count k uses numpy.random.default_rng([SEED, k, d]), the same
scheme as gaucher_downsample.py, so any single draw can be regenerated in isolation.
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
KS = [5, 6]                     # the two configurations quoted in the Table S13 legend
NTOTAL = 31                     # subsample size the bootstrap is used to reach
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

    out = {"seed": SEED, "ndraws": NDRAWS, "n_total": NTOTAL, "route": ROUTE, "nlv": NLV, "k": {}}
    for k in KS:
        nctrl = NTOTAL - k
        t = time.time()
        q2s = np.empty(NDRAWS); aucs = np.empty(NDRAWS); dups = np.empty(NDRAWS)
        for d in range(NDRAWS):
            rng = np.random.default_rng([SEED, k, d])
            sc = rng.choice(cases, size=k, replace=False)
            so = rng.choice(ctrl, size=nctrl, replace=True)      # controls with replacement
            dups[d] = nctrl - len(np.unique(so))
            sel = np.concatenate([sc, so])
            Xs, ys = Xc[sel], y[sel]
            yh = gp.loo_nested(Xs, ys, ROUTE, groups, reps_u, nlv=NLV)
            q2s[d] = gp.q2(ys, yh); aucs[d] = roc_auc_score(ys, yh)
        out["k"][str(k)] = {"n_draws": NDRAWS, "n_controls": nctrl, "q2": q2s.tolist(),
                            "auc": aucs.tolist(), "duplicate_controls": dups.tolist()}
        print("k=%-2d  n=%d (%d controls drawn from 20, %.2f duplicates on average)  "
              "Q2 median %+.4f  AUC median %.4f  Q2>0 %5.1f%%  [%.0f s]"
              % (k, NTOTAL, nctrl, dups.mean(), np.median(q2s), np.median(aucs),
                 100 * np.mean(q2s > 0), time.time() - t), flush=True)

    path = _os.path.join(_RES, "gaucher_bootstrap.json")
    json.dump(out, open(path, "w"))
    print("done after %.0f s -> %s" % (time.time() - t0, path))


if __name__ == "__main__":
    main()
