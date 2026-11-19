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
import sys, pickle, os, numpy as np, nested_generic as ng, altorder as ao, unsup_rep as ur, ordinary_perm as op
from scipy.stats import rankdata
short = sys.argv[1]; sheet = sys.argv[2]
STORE = (_WORK + '/ordinary_perm_results.pkl')
results = pickle.load(open(STORE,'rb')) if os.path.exists(STORE) else {}
if short in results:
    print(f"{short}: already present, skipping."); sys.exit(0)
D0 = ng.build('CS1-4_H.I.'); groups = ao.get_groups(D0['X']); reps_u = ur.centrality_reps(D0['X'], groups)
D = ng.build(sheet); y = D['Y']; N = D['N']; R = rankdata(D['X'], axis=0)
cols1, eff1, pv1, _ = ao.alt_select(groups, R, y, N)
r1 = op.ordinary_perm_test(D['X'], y, cols1, 1, n_perm=1000, seed=201)
ecomb, pcomb = ng.eff_filter(R, y, N); keep = np.where((ecomb>0.3)&(pcomb<=0.05))[0]
rl = ng.group_reps(D['X'][:,keep], R[:,keep], ecomb[keep], N); cols2 = keep[rl]
r2 = op.ordinary_perm_test(D['X'], y, cols2, 1, n_perm=1000, seed=202)
colsv = ng._viptop(D['X'], y, cols1, 15) if len(cols1) > 15 else cols1
rv = op.ordinary_perm_test(D['X'], y, colsv, 1, n_perm=1000, seed=203)
colsu, effu, pvu = ur.unsup_select(reps_u, R, y, N)
ru = op.ordinary_perm_test(D['X'], y, colsu, 1, n_perm=1000, seed=204)
results[short] = dict(A1=r1[:4], A2=r2[:4], VIP=rv[:4], A1u=ru[:4])
pickle.dump(results, open(STORE, 'wb'))
print(f"{short} A1 q={r1[0]:+.3f}(p={r1[2]:.3f}) a={r1[1]:.3f}(p={r1[3]:.3f}) | "
      f"A2 q={r2[0]:+.3f}(p={r2[2]:.3f}) a={r2[1]:.3f}(p={r2[3]:.3f}) | "
      f"VIP q={rv[0]:+.3f}(p={rv[2]:.3f}) a={rv[1]:.3f}(p={rv[3]:.3f}) | "
      f"A1u q={ru[0]:+.3f}(p={ru[2]:.3f}) a={ru[1]:.3f}(p={ru[3]:.3f})")
