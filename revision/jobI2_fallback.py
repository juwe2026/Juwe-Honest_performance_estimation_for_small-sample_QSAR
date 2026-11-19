# CORRECTED SEEDS: the first version of this run used the seed 702 for all routes,
# but 702 is the seed of Route B. The published runs use
# 31 (A1, prun_alt.py), 47 (A1u, prun_uns.py) and 7 (A2, prun.py). The assignment was
# checked with check_seeds.py against the first null values of the deposited files
# (deviation 0). Only with the same seed is the difference of the p-values the effect of the
# correction and not additional Monte Carlo noise.
#
# Own finding: the fallback prediction ytr.mean() in the empty-selection branch reveals the held-out label.
# Here: permutation p-values with the neutral constant 0.5 instead of ytr.mean(), otherwise identical pipeline.
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
import sys,os,json,time
sys.path.insert(0,(_ROOT + '/code'))      # corrected copy first
sys.path.insert(1,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())     # sets sys.path again
for _m in ['nested_generic','altorder','unsup_rep','core']:
    sys.modules.pop(_m,None)
sys.path.insert(0,(_ROOT + '/code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur):
    assert (_ROOT + '/code') in m.__file__, m.__file__
from scipy.stats import rankdata
NP=int(os.environ.get('NP','2000'))
ROUTES=os.environ.get('ROUTES','A1u,A1,A2').split(',')
SEED={'A1':31,'A1u':47,'A2':7}   # as in the published runs
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
out={}
p=os.path.join(_OUT,'I_fallback_seeds.json')
if os.path.exists(p): out=json.load(open(p))
out.setdefault('meta',dict(nperm=NP,seed={'A1':31,'A1u':47,'A2':7},change='empty-selection fallback prediction changed from mean of the training labels to the constant 0.5'))
for k in KEYS:
    y=Y[k]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    out.setdefault(k,{})
    fns={'A1':lambda yy: ao.alt_nested_cv(D1,G,yy,1)[2],
         'A1u':lambda yy: ur.unsup_nested_cv(D1,REPS_U,yy,1)[2],
         'A2':lambda yy: ng.nested_cv(D1,yy,1)[2]}
    for nm in ROUTES:
        if nm in out[k]: print('cache',k,nm,flush=True); continue
        f=fns[nm]; t=time.time()
        yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
        for j in range(NP):
            yp=np.random.default_rng([SEED[nm],j]).permutation(y); h=f(yp)
            qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
        out[k][nm]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),nperm=NP,seconds=round(time.time()-t,1))
        print('%s %-4s Q2 %6.3f p %.4f AUC %.3f p %.4f  %.0fs'%(k,nm,oq,out[k][nm]['p_q2'],oa,out[k][nm]['p_auc'],time.time()-t),flush=True)
        json.dump(out,open(p,'w'),indent=1,default=float)
print('done',ROUTES,flush=True)
