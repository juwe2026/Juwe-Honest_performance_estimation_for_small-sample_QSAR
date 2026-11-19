# -*- coding: utf-8 -*-
"""Morgan fingerprint matrix of the 31 monoterpenoids, read by jobD_rf.py, jobD2_rf_morgan.py, jobE_misc.py,
jobE2_morgan_fixed.py and verify_morgan_fb.py as morgan2_bits.npy.

Morgan fingerprints of radius 2 folded to 2,048 bits (RDKit), computed from the SMILES column of Additional file 3
(sheet Cleaned_Steps1-4, in its row order), of which the 295 bits that vary across the 31 compounds are kept, in bit
order. The file used for the revision was not archived with its generating step. This script rebuilds it, and with
the rebuilt matrix the observed values in morgan2_fixed.json and D_rf_morgan.json are reproduced exactly.
Usage: python revision/make_morgan_bits.py
Writes morgan2_bits.npy to results/revision/ (or QSAR_REV_OUT) and reports whether it equals the file found there.
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
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
XLSX = _osp.environ.get('QSAR_ROUTEA_XLSX', _osp.path.join(_ROOT, 'data', '20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
smiles = pd.read_excel(XLSX, sheet_name='Cleaned_Steps1-4', usecols=['Monoterpenoid', 'SMILES'])['SMILES']
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
F = np.zeros((len(smiles), 2048))
for i, s in enumerate(smiles):
    DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(s)), F[i])
MB = F[:, F.std(axis=0) > 0]
assert MB.shape == (31, 295), MB.shape
out = _osp.path.join(_OUT, 'morgan2_bits.npy')
if _osp.path.exists(out):
    print('equal to the existing file:', bool(np.array_equal(np.load(out), MB)))
_osp.makedirs(_OUT, exist_ok=True); np.save(out, MB); print('written', out, MB.shape)
