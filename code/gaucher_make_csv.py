# -*- coding: utf-8 -*-
"""Rebuild the two Gaucher input files of this archive from the original data of Smit et al. (2007).

Source: gaucherdata.zip of the Biosystems Data Analysis Group, Universiteit van Amsterdam
(https://bdagroup.nl/?p=686), kept unchanged in data/gaucher_source/ together with its licence.
The file Gaucherdata.mat in it holds X (40 spectra x 590 m/z variables), moverz (the 590 m/z
values) and label (20 'control', then 20 'gaucher').

Steps, in this order:
  1. The spectrum in position 38 of the file (patient 18) is left out, which gives the
     39 samples analysed here (20 controls, 19 patients).
  2. median_normalised: every spectrum is divided by its own median over the 590 variables.
  3. autoscaled: every variable of the 39 median-normalised spectra is centred on its mean
     and divided by its standard deviation (n - 1 in the denominator). The analyses of this
     archive read only the median-normalised file and autoscale inside every fold.
  4. Values and m/z column names are written with 10 significant digits, lines end in CR LF.
     Samples are named control_01 to control_20 and gaucher_01 to gaucher_20 by their
     position in the file, so gaucher_18 is the one missing.

Usage:  python code/gaucher_make_csv.py
Writes both files to results/recalc/ (or QSAR_RECALC_DIR) and compares them byte for byte with
data/Gaucherdata_paper39_median_normalised.csv and data/Gaucherdata_paper39_autoscaled.csv.
The zip file can be given elsewhere with QSAR_GAUCHER_ZIP.
"""
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)

import hashlib, io, zipfile
import numpy as np
import scipy.io as sio

ZIP = _osp.environ.get('QSAR_GAUCHER_ZIP', _osp.path.join(_ROOT, 'data', 'gaucher_source', 'gaucherdata.zip'))
DROP = 37          # zero-based position of patient 18 in Gaucherdata.mat

with zipfile.ZipFile(ZIP) as z:
    m = sio.loadmat(io.BytesIO(z.read('Gaucherdata.mat')))
X = np.asarray(m['X'], float)
mz = np.asarray(m['moverz'], float).ravel()
lab = [str(np.asarray(v).ravel()[0]) for v in m['label'].ravel()]
assert X.shape == (40, 590) and mz.size == 590 and lab == ['control'] * 20 + ['gaucher'] * 20

names = ['control_%02d' % (i + 1) for i in range(20)] + ['gaucher_%02d' % (i + 1) for i in range(20)]
keep = [i for i in range(40) if i != DROP]
Xm = (X / np.median(X, axis=1, keepdims=True))[keep]
Z = (Xm - Xm.mean(axis=0)) / Xm.std(axis=0, ddof=1)

def write(fn, A):
    rows = ['sample,class,' + ','.join('%.10g' % v for v in mz)]
    for j, i in enumerate(keep):
        rows.append(names[i] + ',' + lab[i] + ',' + ','.join('%.10g' % v for v in A[j]))
    with open(fn, 'w', newline='') as f:
        f.write('\r\n'.join(rows) + '\r\n')
    return fn

for tag, A in (('median_normalised', Xm), ('autoscaled', Z)):
    new = write(_osp.path.join(_WORK, 'Gaucherdata_paper39_%s.csv' % tag), A)
    ref = _osp.path.join(_ROOT, 'data', 'Gaucherdata_paper39_%s.csv' % tag)
    h = [hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in (new, ref)]
    print('%-18s %s  %s' % (tag, h[0][:16], 'byte-identical with data/' if h[0] == h[1] else 'DIFFERS from data/'))
