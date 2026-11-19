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
import sys, json, time, numpy as np, pandas as pd
from scipy.stats import spearmanr
def gruppieren_orig(A):          # literally identical with datasets.gruppieren (benchmark)
    P = np.corrcoef(A, rowvar=False); S = spearmanr(A).statistic
    R = (np.abs(np.nan_to_num(P)) > 0.7) | (np.abs(np.nan_to_num(S)) > 0.7); np.fill_diagonal(R, False)
    var = A.var(axis=0, ddof=0); offen = set(range(A.shape[1])); gruppen = []
    while offen:
        grade = {i: sum(1 for j in offen if j != i and R[i, j]) for i in offen}
        mx = max(grade.values()); kand = [i for i in offen if grade[i] == mx]
        hub = max(kand, key=lambda i: (var[i], -i))
        mit = {j for j in offen if j != hub and R[hub, j]}
        gruppen.append([hub] + sorted(mit)); offen -= ({hub} | mit)
    return gruppen
def gruppieren_vec(A):           # vectorised, same rule: degree among the open ones, then largest variance, then smallest index
    P = np.corrcoef(A, rowvar=False); S = spearmanr(A).statistic
    R = (np.abs(np.nan_to_num(P)) > 0.7) | (np.abs(np.nan_to_num(S)) > 0.7); np.fill_diagonal(R, False)
    var = A.var(axis=0, ddof=0); k = A.shape[1]; offen = np.ones(k, bool); gruppen = []
    deg = R.sum(1).astype(np.int64)
    while offen.any():
        d = np.where(offen, deg, -1); mx = d.max(); kand = np.where(d == mx)[0]
        vm = var[kand].max(); kand = kand[var[kand] == vm]; hub = int(kand.min())
        mit = np.where(R[hub] & offen)[0]; mit = mit[mit != hub]
        gruppen.append([hub] + sorted(int(j) for j in mit))
        weg = np.r_[hub, mit]; offen[weg] = False
        deg -= R[:, weg].sum(1)
    return gruppen
if __name__ == '__main__':
    X = pd.read_csv((_ROOT + '/revision/lone/desc_pool.csv')); X = X[[c for c in X.columns if c != 'No']].values.astype(float)
    t = time.time(); g2 = gruppieren_vec(X); t2 = time.time() - t
    t = time.time(); g1 = gruppieren_orig(X); t1 = time.time() - t
    print('peptides orig', len(g1), f'{t1:.1f}s', 'vec', len(g2), f'{t2:.2f}s', 'identical', g1 == g2, flush=True)
    for ds in ['DS01', 'DS07', 'DS16']:
        Xd = pd.read_csv(f'{_ROOT}/benchmark/datasets_k/{ds}_X.csv'); cols = list(Xd.columns)
        info = json.load(open(f'{_ROOT}/benchmark/datasets_k/{ds}_info.json'))
        ref = [[cols.index(c) for c in grp] for grp in info['gruppen']]
        gv = gruppieren_vec(Xd.values.astype(float))
        print(ds, 'saved', len(ref), 'vec', len(gv), 'identical', ref == gv, flush=True)
