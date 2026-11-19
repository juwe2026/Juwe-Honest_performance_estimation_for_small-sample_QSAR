# Route A2, H. influenzae: separate permutation run in the deposited state.
# The four other strains come from prun.py (seed 7), H. influenzae from
# perm_runner.py (seed 42, via nested_h.py). Checked with check_seed_a2c.py against the
# first null value of perm_HI.npz: seed 42 matches exactly, 7 does not.
# Here recomputed with corrected fallback branch and 2,000 permutations.
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
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX',(_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,(_ROOT + '/code'))
import nested_generic as ng
assert (_ROOT + '/code') in ng.__file__
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
NP=2000; SEED=42
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
y=Y['HI']; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
f=lambda yy: ng.nested_cv(D1,yy,1)[2]
t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
for j in range(NP):
    yp=np.random.default_rng([SEED,j]).permutation(y); h=f(yp)
    qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
qs=np.array(qs); as_=np.array(as_)
res=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),
         nperm=NP,seed=SEED,seconds=round(time.time()-t,1))
json.dump(res,open((_OUT + '/A2_HI_seed42.json'),'w'),indent=1,default=float)
print('A2 HI seed 42: Q2 %.4f p %.4f | AUC %.4f p %.4f  (%.0fs)'%(oq,res['p_q2'],oa,res['p_auc'],res['seconds']),flush=True)
