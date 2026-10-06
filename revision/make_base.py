# -*- coding: utf-8 -*-
"""Writes base.pkl, the shared input of the revision scripts (loaded by common.py).

base.pkl holds a dict with
  D     ng.build('Cleaned_Steps1-4', ycol=...) for the five strains (keys HI, SA, SPneu, SPyo, PA)
  g     the 182 correlation groups of altorder.build_groups
  ADJ   the adjacency matrix of that grouping (|r| or |rho| > 0.7)
  reps  the label-free representatives of unsup_rep.centrality_reps (Route A1u)
built from the archive code and the Route A matrix (Additional file 6, placed as described in the
README). It is not bundled because of its size (about 175 MB) and because this script rebuilds it.

Usage:  python revision/make_base.py [OUT]       default OUT: results/recalc/base.pkl
common.py expects the file at /home/claude/rev/base.pkl (see README_revision.md on paths).
"""
import os, sys, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur

COLS = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPneu': 'S. Pneu (S. pneumoniae)',
        'SPyo': 'S. pyo (S. pyogenes)', 'PA': 'P. A. (P. aeruginosa)'}
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ng._RES, 'base.pkl')

D = {k: ng.build('Cleaned_Steps1-4', ycol=c) for k, c in COLS.items()}
X = D['HI']['X']
g, ADJ = ao.build_groups(X)
assert len(g) == 182, len(g)
reps = ur.centrality_reps(X, g)
pickle.dump(dict(D=D, g=g, ADJ=ADJ, reps=reps), open(OUT, 'wb'))
print('written', OUT, '-', X.shape[1], 'descriptors,', len(g), 'groups')
