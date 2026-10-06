# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
"""Route B, 12-descriptor set (C_Count removed, perfectly collinear with
EsterLacton_flag): fully leave-one-out nested permutation test, one and three
latent variables, 2,000 permutations per strain per model.

This fills the one archive gap flagged in the independent review: the raw
permutation outputs for every route except the published 12-descriptor Route B
model were archived, while that model itself was represented only by a cached
summary object (figures/data_recalc/routeB12.json). Running this script
reproduces that summary from scratch (the observed Q2/AUC values it prints
match routeB12.json to the last printed digit) and, unlike the earlier
computation, also saves every individual permutation's Q2 and AUC to
results/raw_permutations/, in the same format dataio.perm() already reads for
every other route (permB_12_<strain>.npz and permB_12nlv3_<strain>.npz).

Same seed convention as code/run_routeB.py (base seed 702 for the reduced
descriptor set): permutation k for strain s draws
np.random.default_rng([702, k]).permutation(y). Both models (nLV 1 and nLV 3)
are evaluated on the SAME permuted label vector for a given k, exactly as
"2000 Permutationen je Modell, Seed 702" in routeB12.json's own note already
implied, so the permutation identity lines up across nLV 1 and nLV 3 the same
way it lines up with the other 25-/13-descriptor permutation files.

Usage: python run_routeB12.py [n_permutations] [outdir]
  n_permutations defaults to 2000; outdir defaults to ../results/raw_permutations
"""
import os, sys, time
import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
RAW_DEFAULT = os.path.join(REPO_ROOT, "results", "raw_permutations")

N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else RAW_DEFAULT

SHORT = {"H.I.": "HI", "S.A.": "SA", "S.Pneu": "SPneu", "S.Pyo": "SPyo", "P.A.": "PA"}
YCOL = {"H.I.": "H. I. (H. influenzae)", "S.A.": "S. A. (Staph. Aureus)",
        "S.Pneu": "S. Pneu (S. pneumoniae)", "S.Pyo": "S. pyo (S. pyogenes)",
        "P.A.": "P. A. (P. aeruginosa)"}


def loo(X, y, nlv):
    N = len(y)
    yh = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xt = X[tr]; mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd == 0] = 1
        m = PLSRegression(n_components=min(nlv, X.shape[1], N - 2), scale=False).fit((Xt - mu) / sd, y[tr])
        yh[i] = m.predict((X[i:i + 1] - mu) / sd).ravel()[0]
    return yh


def q2(y, yh):
    return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)


def main():
    full = pd.read_pickle(os.path.join(RAW_DEFAULT, "desc_full.pkl"))
    reduced13 = list(np.load(os.path.join(RAW_DEFAULT, "reduced_desc.npy"), allow_pickle=True))
    cols12 = [c for c in reduced13 if c != "C_Count"]
    assert len(cols12) == 12, cols12
    os.makedirs(OUTDIR, exist_ok=True)

    for s, short in SHORT.items():
        t0 = time.time()
        y = (full[YCOL[s]].astype(str).str.strip().str.lower() == "active").astype(float).values
        X = full[cols12].astype(float).values

        obs = {}
        for nlv in (1, 3):
            yh = loo(X, y, nlv)
            obs[nlv] = (q2(y, yh), roc_auc_score(y, yh))

        q2s = {1: [], 3: []}
        aucs = {1: [], 3: []}
        for k in range(N_PERM):
            yp = np.random.default_rng([702, k]).permutation(y)
            for nlv in (1, 3):
                yhp = loo(X, yp, nlv)
                q2s[nlv].append(q2(yp, yhp))
                aucs[nlv].append(roc_auc_score(yp, yhp))
            if (k + 1) % 200 == 0:
                print("  %s: %d/%d permutations, %.0fs" % (s, k + 1, N_PERM, time.time() - t0), flush=True)

        for nlv, tag in ((1, "permB_12"), (3, "permB_12nlv3")):
            oq, oa = obs[nlv]
            qs = np.array(q2s[nlv]); as_ = np.array(aucs[nlv])
            n = len(qs)
            pq = (1 + np.sum(qs >= oq - 1e-9)) / (1 + n)
            pa = (1 + np.sum(as_ >= oa - 1e-9)) / (1 + n)
            np.savez(os.path.join(OUTDIR, "%s_%s.npz" % (tag, short)),
                     q2s=qs, aucs=as_, obs_q2=oq, obs_auc=oa)
            print("%s %s nLV=%d: n=%d | Q2=%.4f AUC=%.4f | p(Q2)=%.4f p(AUC)=%.4f | n_ge_q2=%d n_ge_auc=%d"
                  % (s, tag, nlv, n, oq, oa, pq, pa, int(np.sum(qs >= oq - 1e-9)), int(np.sum(as_ >= oa - 1e-9))))
        print("  %s done in %.0fs" % (s, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
