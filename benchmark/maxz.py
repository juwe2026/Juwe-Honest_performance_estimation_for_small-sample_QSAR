"""Label-free criterion R2 for numerical instability: for every held-out compound (LOO) the largest |z| over all
cleaned descriptors, autoscaled with mean and SD (ddof=0) of the training compounds. Descriptors
that are constant in the training compounds (SD < 1e-6) are not evaluated. No label is used.
Reproduces ergebnisse/maxz_labelfrei.json of the analysis under the original rule exactly (checked 25.10.2026).
Checked again on 6 November 2026: ergebnisse_k/ and ergebnisse_m/ byte for byte, ergebnisse/ to within a relative
difference of 4e-13 (no series changes its R2 flag).
Usage: python3 maxz.py <datasets folder> <output file>
"""
import json, sys, glob
import numpy as np, pandas as pd

def maxz(A):
    out = []
    for i in range(A.shape[0]):
        tr = np.delete(A, i, 0)
        mu = tr.mean(0); sd = tr.std(0, ddof=0); ok = sd >= 1e-6
        out.append(float(np.max(np.abs((A[i, ok] - mu[ok]) / sd[ok]))))
    return out

if __name__ == '__main__':
    ordner, aus = sys.argv[1], sys.argv[2]
    dss = sorted({f.split('/')[-1][:4] for f in glob.glob(f'{ordner}/DS*_X.csv')})
    res = {ds: maxz(pd.read_csv(f'{ordner}/{ds}_X.csv').values.astype(float)) for ds in dss}
    json.dump(res, open(aus, 'w'))
