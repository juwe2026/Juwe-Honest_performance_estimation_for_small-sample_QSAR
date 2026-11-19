# Which seed lies behind which deposited permutation file?
# Check: recompute the first permutations with the ORIGINAL CODE (faulty fallback branch)
# and the presumed seed and compare them with the deposited null values.
# If they agree, the seed assignment is right, and the corrected runs must use the same
# seed, otherwise Monte Carlo noise enters the difference.
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
import sys, os, numpy as np
os.environ.setdefault('QSAR_ROUTEA_XLSX',(_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
import nested_generic as ng, altorder as ao, unsup_rep as ur
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
assert _osp.path.realpath(_osp.path.dirname(ng.__file__)) == _osp.path.realpath(_ROOT + '/code'), ng.__file__

FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))

RP=(_ROOT + '/results/raw_permutations/')
ARMS={
 'perm_alt_%s.npz'  : (31, lambda D1,G,yy: ao.alt_nested_cv(D1,G,yy,1)[2]),
 'perm_uns_%s.npz'  : (47, None),
 'perm_%s.npz'      : (7,  lambda D1,G,yy: ng.nested_cv_a2(D1,yy,1)[2] if hasattr(ng,'nested_cv_a2') else None),
 'perm_a1e_%s.npz'  : (59, lambda D1,G,yy: ao.alt_nested_cv(D1,G,yy,1,15)[2]),
 'perm_a1vip1_%s.npz':(137,lambda D1,G,yy: ao.alt_viptop(D1,G,yy,1,15,1)[2]),
 'perm15_%s.npz'    : (11, lambda D1,G,yy: ng.nested_cv_a2_topk(D1,yy,1,15)[2]),
}
K='HI'; y=Y[K]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
for pat,(seed,fn) in ARMS.items():
    if fn is None: continue
    f=RP+pat%K
    if not os.path.exists(f):
        print('%-22s file missing'%pat); continue
    z=np.load(f,allow_pickle=True); q2s=np.asarray(z['q2s'],float)
    ok=[]
    for j in range(3):
        yp=np.random.default_rng([seed,j]).permutation(y)
        h=fn(D1,G,yp)
        if h is None: ok=None; break
        ok.append(abs(q2(yp,h)-q2s[j]))
    if ok is None:
        print('%-22s function not available'%pat); continue
    print('%-22s seed %4d | deviation of the first three null values: %s'
          %(pat,seed,', '.join('%.2e'%d for d in ok)))
