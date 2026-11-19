# -*- coding: utf-8 -*-
"""Do the nested Gaucher results depend on which spectrum is left out?

Smit et al. (2007, section 4.1) left out one patient spectrum that appeared flawed, without saying which one in the
published file. The analyses of this archive leave out spectrum 38 of Gaucherdata.mat (patient 18, see
gaucher_spectrum_check.py). This script repeats the nested leave-one-out analysis of Routes A1, A1u and A2 with one
latent variable (../../code/gaucher_pipe.py, as for Table S21) on
  - the archived input data/Gaucherdata_paper39_median_normalised.csv,
  - all 40 spectra,
  - the 39 spectra left after leaving out each one of the 20 patient spectra in turn.
Every variant is median-normalised as in Smit et al., and the correlation groups and label-free representatives are
rebuilt on its own spectra. Only observed values are computed, no permutations.
Usage: python gaucher_selection_sensitivity.py <archive_root> <out.json>
"""
import io, json, os, sys, zipfile
import numpy as np, scipy.io as sio
from sklearn.metrics import roc_auc_score
ROOT, OUT = sys.argv[1:3]
sys.path.insert(0, os.path.join(ROOT, 'code')); import gaucher_pipe as gp
m = sio.loadmat(io.BytesIO(zipfile.ZipFile(os.path.join(ROOT, 'data', 'gaucher_source', 'gaucherdata.zip')).read('Gaucherdata.mat')))
X40 = np.asarray(m['X'], float); X40 = X40 / np.median(X40, axis=1, keepdims=True)
mz40 = np.asarray(m['moverz'], float).ravel(); y40 = np.r_[np.zeros(20), np.ones(20)]
res = {}
for v in ['archived_input', 'all_40'] + ['without_gaucher_%02d' % k for k in range(1, 21)]:
    if v == 'archived_input':
        X, y, mz, _ = gp.load(os.path.join(ROOT, 'data', 'Gaucherdata_paper39_median_normalised.csv'))
    else:
        keep = list(range(40)) if v == 'all_40' else [i for i in range(40) if i != 19 + int(v[-2:])]
        X, y, mz = X40[keep], y40[keep], mz40
    Xc, _, _ = gp.clean(X, mz); groups, _ = gp.build_groups(Xc); reps_u = gp.centrality_reps(Xc, groups)
    r = dict(n=int(len(y)), patients=int(y.sum()), groups=len(groups))
    for route in ('A1', 'A1u', 'A2'):
        yh = gp.loo_nested(Xc, y, route, groups, reps_u, nlv=1)
        r[route] = dict(q2=float(gp.q2(y, yh)), auc=float(roc_auc_score(y, yh)))
    res[v] = r
    print('%-20s n=%d groups=%d  ' % (v, r['n'], r['groups']) +
          '  '.join('%s Q2 %.3f AUC %.3f' % (k, r[k]['q2'], r[k]['auc']) for k in ('A1', 'A1u', 'A2')), flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
