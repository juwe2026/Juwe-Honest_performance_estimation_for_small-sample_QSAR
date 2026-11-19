# Permutation p-value with treatment of ties (revision of 03.10.2026).
# p = (1 + #{null values >= observed value}) / (1 + number of permutations).
# AUC values are fractions k/(n1*n0). Equal fractions are stored by the computer, depending on
# the way of computation, with deviations of about 1e-16. The former strict comparison
# did not always count such ties and made p slightly too small.
# Now a null value counts as >= observed if it lies at most 1e-9 below it.
# Every call stores the null distribution under out/nulls/ for checking.
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
import os as _os, sys as _sys, numpy as _np
_TIE_EPS = 1e-9
_NULLDIR = (_OUT + '/nulls')
_CNT = [0]
def pval(null, obs, tag=None):
    null = _np.asarray(null, float); obs = float(obs)
    strict = float((1 + _np.sum(null >= obs)) / (1 + len(null)))
    tie = float((1 + _np.sum(null >= obs - _TIE_EPS)) / (1 + len(null)))
    try:
        _os.makedirs(_NULLDIR, exist_ok=True); _CNT[0] += 1
        nm = '%s_%s_%d_%04d' % (_os.path.basename(_sys.argv[0]).replace('.py', ''),
                               '-'.join(_sys.argv[1:]) or 'x', _os.getpid(), _CNT[0])
        _np.savez_compressed(_os.path.join(_NULLDIR, nm + '.npz'), null=null, obs=obs, strict=strict, tie=tie,
                             tag=str(tag), argv=' '.join(_sys.argv))
    except Exception:
        pass
    return tie
