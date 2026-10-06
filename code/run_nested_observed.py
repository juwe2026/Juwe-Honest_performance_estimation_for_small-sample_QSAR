"""Observed fully nested leave-one-out Q2 and AUC of Routes A1, A1u and A2 for all five strains.

Added for the revision (October 2026) as the entry point of the nested validation. `nested_generic.py`
holds the routines but has no command-line entry point of its own. This script calls the same
functions that the whole-pipeline permutation scripts use for their observed values:

  Route A1   altorder.alt_nested_cv     (correlation grouping of the whole pool first; representative
                                         chosen by effect size inside every fold, then the filter)
  Route A1u  unsup_rep.unsup_nested_cv  (representative chosen by centrality, without labels)
  Route A2   nested_generic.nested_cv   (filter on all 2,208 descriptors first, then grouping of the
                                         survivors inside every fold)

One latent variable throughout, as in the article (Table 1, blocks a and b). Needs the Route A
descriptor matrix in data/ (see README, "One input has to be placed by hand").

Usage:  python code/run_nested_observed.py
Writes results/recalc/nested_observed.json and prints one line per strain and route.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nested_generic as ng, altorder as ao, unsup_rep as ur

YCOL = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPneu': 'S. Pneu (S. pneumoniae)',
        'SPyo': 'S. pyo (S. pyogenes)', 'PA': 'P. A. (P. aeruginosa)'}
OUT = os.path.join(ng._RES, 'nested_observed.json')

res = {}
for short, col in YCOL.items():
    D = ng.build('Cleaned_Steps1-4', ycol=col)
    groups = ao.get_groups(D['X'])
    reps = ur.centrality_reps(D['X'], groups)
    q1, a1, _ = ao.alt_nested_cv(D, groups, D['Y'], 1)
    qu, au, _ = ur.unsup_nested_cv(D, reps, D['Y'], 1)
    q2, a2, _ = ng.nested_cv(D, D['Y'], 1)
    res[short] = {'A1': {'q2': q1, 'auc': a1}, 'A1u': {'q2': qu, 'auc': au}, 'A2': {'q2': q2, 'auc': a2},
                  'n_groups': len(groups), 'n_active': int(D['Y'].sum())}
    for r in ('A1', 'A1u', 'A2'):
        print('%-6s %-4s nested Q2 = %6.3f   AUC = %.3f' % (short, r, res[short][r]['q2'], res[short][r]['auc']))
json.dump(res, open(OUT, 'w'), indent=1)
print('written', OUT)
