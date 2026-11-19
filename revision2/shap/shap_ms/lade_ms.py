"""Data sets of the manuscript (Revision 1, Additional files 4 and 6 of Revision 1, state of 06.10.2026)."""
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
# QSAR_REV1_AF_DIR: folder with Additional_file_6_RouteA_descriptor_matrix 06.10.2026.xlsx and
# Additional_file_4_second_MIC_data_set 06.10.2026.xlsx of the first revision (not part of the archive).
P = _osp.environ.get('QSAR_REV1_AF_DIR', _WORK) + '/'
def monoterpene():
    df = pd.read_excel(P + 'Additional_file_6_RouteA_descriptor_matrix 06.10.2026.xlsx', sheet_name='Cleaned_Steps1-4')
    ep = {'HI': 'H. I. (H. influenzae)', 'PA': 'P. A. (P. aeruginosa)', 'SA': 'S. A. (S. aureus)',
          'SPneu': 'S. Pneu (S. pneumoniae)', 'SPyo': 'S. pyo (S. pyogenes)'}
    X = df.iloc[:, 7:].astype(float)
    Y = {k: (df[v].str.strip() == 'active').astype(float).values for k, v in ep.items()}
    return X, Y
def peptide():
    d = pd.read_excel(P + 'Additional_file_4_second_MIC_data_set 06.10.2026.xlsx', sheet_name='Descriptors_1275')
    c = pd.read_excel(P + 'Additional_file_4_second_MIC_data_set 06.10.2026.xlsx', sheet_name='Peptides_MIC_classes')
    assert list(d['No'].astype(str)) == list(c['No'].astype(str))
    X = d.drop(columns=['No']).astype(float)
    Y = {k: (c[f'Class rule A {k}'].str.strip() == 'active').astype(float).values for k in ['EC', 'SA', 'PA', 'AB']}
    return X, Y
if __name__ == '__main__':
    for f in (monoterpene, peptide):
        X, Y = f(); print(f.__name__, X.shape, {k: int(v.sum()) for k, v in Y.items()}, 'NaN', int(X.isna().sum().sum()))
