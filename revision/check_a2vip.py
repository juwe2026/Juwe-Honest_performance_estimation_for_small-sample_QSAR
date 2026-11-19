# Check of the reconstruction: does nested_cv_a2_viptop reproduce the deposited
# observed values of perm_vip_*.npz? The fallback branch is corrected in fixcode,
# but it practically never applies with the true labels, so the
# observed values must agree if the reconstruction is correct.
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
import sys, os, glob
os.environ.setdefault('QSAR_ROUTEA_XLSX',(_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,(_ROOT + '/code'))
import nested_generic as ng, a2vip
from scipy.stats import rankdata
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
RP=(_ROOT + '/results/raw_permutations/')
for nlv_vip in (2,1):
    print('--- nlv_vip =',nlv_vip)
    for k in KEYS:
        y=Y[k]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
        q,a,_=a2vip.nested_cv_a2_viptop(D1,y,1,15,nlv_vip)
        z=np.load(RP+'perm_vip_%s.npz'%k,allow_pickle=True)
        oq=float(z['obs_q2']); oa=float(z['obs_auc'])
        print('%-6s Q2 %7.4f (archive %7.4f, d=%7.4f) | AUC %6.4f (archive %6.4f, d=%6.4f)'
              %(k,q,oq,q-oq,a,oa,a-oa))
