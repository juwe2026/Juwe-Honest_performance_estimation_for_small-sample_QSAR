# Revision 03.10.2026: permutation p-values count ties (tolerance 1e-9), see CHANGELOG_revision.md
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
import sys, os, numpy as np, nested_h as nh

STORE = (_WORK + '/perm_HI.npz')
NLV = 1
batch = int(sys.argv[1]) if len(sys.argv) > 1 else 400

# store the observed value once
if os.path.exists(STORE):
    d = np.load(STORE)
    q2s = list(d['q2s']); aucs = list(d['aucs'])
    obs_q2 = float(d['obs_q2']); obs_auc = float(d['obs_auc'])
else:
    obs_q2, obs_auc, _, nf = nh.nested_cv(nh.Y, nlv=NLV)
    q2s, aucs = [], []
    print(f'Observed: Q2={obs_q2:.4f} AUC={obs_auc:.4f}')

start = len(q2s)
for k in range(start, start + batch):
    rng = np.random.default_rng([42, k])          # reproducible, independent per index
    yv = rng.permutation(nh.Y)
    q, a, _, _ = nh.nested_cv(yv, nlv=NLV)
    q2s.append(q); aucs.append(a)

np.savez(STORE, q2s=np.array(q2s), aucs=np.array(aucs),
         obs_q2=obs_q2, obs_auc=obs_auc)
n = len(q2s)
pQ = (1 + np.sum(np.array(q2s) >= obs_q2 - 1e-9)) / (1 + n)
pA = (1 + np.sum(np.array(aucs) >= obs_auc - 1e-9)) / (1 + n)
print(f'perms total: {n} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | '
      f'perm-mean Q2={np.mean(q2s):.3f} AUC={np.mean(aucs):.3f}')
