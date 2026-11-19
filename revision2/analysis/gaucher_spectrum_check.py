# -*- coding: utf-8 -*-
"""Which of the 40 spectra in Gaucherdata.mat is the one Smit et al. (2007, section 4.1) excluded as flawed?

Smit et al. excluded one patient spectrum on visual inspection but do not say which one in the published file. This
script compares the 40 median-normalised spectra by label-free criteria, each computed for one spectrum against the
other 39: the Spearman correlation with their median spectrum, Hotelling's T2 and the residual Q in a principal
component analysis of the other 39 (autoscaled, five components), and the largest intensity before normalisation.
It also counts the correlation groups of the analysis (gaucher_pipe.build_groups) with all 40 spectra and with each
patient spectrum left out.
Usage: python gaucher_spectrum_check.py <archive_root> <out.json>
"""
import io, json, os, sys, zipfile
import numpy as np, scipy.io as sio
from scipy.stats import spearmanr
ROOT, OUT = sys.argv[1:3]
sys.path.insert(0, os.path.join(ROOT, 'code')); import gaucher_pipe as gp
m = sio.loadmat(io.BytesIO(zipfile.ZipFile(os.path.join(ROOT, 'data', 'gaucher_source', 'gaucherdata.zip')).read('Gaucherdata.mat')))
X = np.asarray(m['X'], float); mz = np.asarray(m['moverz'], float).ravel()
names = ['control_%02d' % (i + 1) for i in range(20)] + ['gaucher_%02d' % (i + 1) for i in range(20)]
Xm = X / np.median(X, axis=1, keepdims=True)
rho, T2, Q = np.zeros(40), np.zeros(40), np.zeros(40)
for i in range(40):
    rest = np.delete(Xm, i, 0)
    rho[i] = spearmanr(Xm[i], np.median(rest, 0))[0]
    mu = rest.mean(0); sd = rest.std(0, ddof=1); sd[sd == 0] = 1
    A = (rest - mu) / sd; _, S, Vt = np.linalg.svd(A, full_matrices=False); P = Vt[:5].T
    lam = S[:5] ** 2 / (len(A) - 1); z = (Xm[i] - mu) / sd; t = z @ P
    T2[i] = float(np.sum(t ** 2 / lam)); Q[i] = float(np.sum((z - t @ P.T) ** 2))
rank = lambda v: (np.argsort(np.argsort(-v)) + 1)
per = {names[i]: dict(spearman_with_median_of_others=round(float(rho[i]), 4), rank_lowest_correlation=int(rank(-rho)[i]),
                      hotelling_T2=round(float(T2[i]), 2), rank_T2=int(rank(T2)[i]), residual_Q=round(float(Q[i]), 1), rank_Q=int(rank(Q)[i]),
                      max_intensity=round(float(X[i].max()), 2), rank_lowest_max=int(rank(-X.max(1))[i]))
       for i in range(40)}
groups = {}
for v in ['all_40'] + ['without_gaucher_%02d' % k for k in range(1, 21)]:
    keep = list(range(40)) if v == 'all_40' else [i for i in range(40) if i != 19 + int(v[-2:])]
    Xc, _, _ = gp.clean(Xm[keep], mz); g, _ = gp.build_groups(Xc); groups[v] = len(g)
res = dict(per_spectrum=per, correlation_groups=groups)
json.dump(res, open(OUT, 'w'), indent=1)
order = sorted(per, key=lambda n: per[n]['rank_lowest_correlation'] + per[n]['rank_T2'] + per[n]['rank_Q'])
for n in order[:5]: print(n, per[n])
print('groups', groups)
