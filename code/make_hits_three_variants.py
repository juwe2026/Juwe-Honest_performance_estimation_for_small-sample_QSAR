# -*- coding: utf-8 -*-
"""Regenerates results/raw_permutations/hits_three_variants.npy (Table S6, Fig. 7, Fig. S6).

Counts how many descriptors pass the univariate filter (effect size > 0.3 and p <= 0.05) under
three reduction orders, once on the observed classes and 1,000 times under permuted classes:

  A2   filter on all 2,208 descriptors
  A1   grouping first, representative of each group chosen by effect size (label-dependent);
       stored under the key "A1s" in the archived object
  A1u  grouping first, representative of each group chosen by centrality (label-free)

The generating script was missing from the submitted archive and was reconstructed for the
revision. Permutation scheme of the archived object: a new generator numpy.random.default_rng(9)
for each strain, and its 1,000 successive calls of permutation(y). With this scheme the script
reproduces all five observed triples and all fifteen null arrays of the archived object element
for element, and the file it writes is byte-identical to the archived one. An earlier draft of this script used default_rng([8201, j]); that scheme does not
reproduce the archived null arrays and was never used for a published number.

Inputs are built from the archive code itself: the Route A matrix (Additional file 6, placed as
described in the README), the 182 correlation groups of altorder.build_groups and the label-free
representatives of unsup_rep.centrality_reps. No intermediate pickle is needed.

Usage:  python code/make_hits_three_variants.py
Writes results/recalc/hits_three_variants.npy (QSAR_RECALC_DIR overrides the folder), compares it
with the archived results/raw_permutations/hits_three_variants.npy and exits non-zero on any
difference.
"""
import os, sys
import numpy as np
from scipy.stats import rankdata
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nested_generic as ng, altorder as ao, unsup_rep as ur

NPERM = 1000
SEED = 9
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCHIVED = os.path.join(ROOT, 'results', 'raw_permutations', 'hits_three_variants.npy')
OUT = os.path.join(ng._RES, 'hits_three_variants.npy')
YCOL = {'H.I.': 'H. I. (H. influenzae)', 'S.A.': 'S. A. (S. aureus)', 'S.Pneu': 'S. Pneu (S. pneumoniae)',
        'S.Pyo': 'S. pyo (S. pyogenes)', 'P.A.': 'P. A. (P. aeruginosa)'}


def hit_counts(R, y, N, groups, reps_unsup):
    """Hit counts of the three variants for one class vector."""
    eff, pv = ng.eff_filter(R, y, N)
    passed = (eff > 0.3) & (pv <= 0.05)
    reps_eff = np.array([g[int(np.argmax(eff[g]))] for g in groups])
    return (int(passed.sum()),
            int(passed[reps_eff].sum()),
            int(passed[reps_unsup].sum()))


def main():
    D = {lab: ng.build('Cleaned_Steps1-4', ycol=col) for lab, col in YCOL.items()}
    X = D['H.I.']['X']
    groups, _ = ao.build_groups(X)
    assert len(groups) == 182, len(groups)
    reps_unsup = np.asarray(ur.centrality_reps(X, groups))
    out = {}
    for lab in YCOL:
        assert np.array_equal(D[lab]['X'], X)          # one descriptor matrix, five class vectors
        y = D[lab]['Y']; N = D[lab]['N']
        R = rankdata(X, axis=0)
        o2, o1, o1u = hit_counts(R, y, N, groups, reps_unsup)
        rng = np.random.default_rng(SEED)              # new generator for every strain
        null = np.array([hit_counts(R, rng.permutation(y), N, groups, reps_unsup) for _ in range(NPERM)],
                        dtype=np.int64)
        out[lab] = dict(obs_A2=o2, null_A2=null[:, 0].copy(), obs_A1s=o1, null_A1s=null[:, 1].copy(),
                        obs_A1u=o1u, null_A1u=null[:, 2].copy())
        print(lab, 'observed A2/A1/A1u:', o2, o1, o1u, flush=True)
    np.save(OUT, out, allow_pickle=True)
    print('written', OUT)
    if os.path.exists(ARCHIVED):
        A = np.load(ARCHIVED, allow_pickle=True).item()
        bad = [(lab, k) for lab in YCOL for k in out[lab]
               if not np.array_equal(np.asarray(out[lab][k]), np.asarray(A[lab][k]))]
        assert set(A) == set(out) and all(set(A[lab]) == set(out[lab]) for lab in YCOL), 'keys differ'
        print('compared with the archived object: %s' % ('identical, all 30 entries' if not bad else 'DIFFERENT %s' % bad))
        sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
