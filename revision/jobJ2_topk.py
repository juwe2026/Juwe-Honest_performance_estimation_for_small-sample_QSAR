# CORRECTED SEEDS: the first version used 702 for all arms (that is the seed
# of Route B). The published runs use 59 (A1 top-15 by effect size,
# prun_a1.py), 137 (A1 top-15 by VIP with one component, prun_vip1.py), 11 (A2
# top-15 by effect size, prun15.py) and 23 (A2 top-15 by VIP, prun_vip.py).
# Checked with check_seeds.py against the first null values of the deposited files
# (deviation 0 for 31, 59, 137 and 11).
#
# Addendum to the code correction: the n:p-constrained selections (Table 3 of Revision 1) and the
# VIP selection (Table 4 of Revision 1) run through the same faulty fallback branch as the
# Routes A1/A1u/A2 and were overlooked in the first correction.
# Here: recomputed with the constant 0.5 AND with the full 2,000 permutations
# (the deposited files perm15_* and perm_vip_* have only 1,900 to 1,982).
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
import nested_generic as ng, altorder as ao, unsup_rep as ur, a2vip
for m in (ng,ao,ur):
    assert (_ROOT + '/code') in m.__file__, m.__file__
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

NP=int(os.environ.get('NP','2000'))
ARMS=os.environ.get('ARMS','A1_top15,A1_viptop1,A1_top6,A2_top15,A2_viptop').split(',')
CACHE=(_OUT + '/cache_topk_seeds'); os.makedirs(CACHE,exist_ok=True)
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
SEED={'A1_top15':59,'A1_viptop1':137,'A2_top15':11,'A2_viptop':23,'A1_viptop':73,'A1_top6':91}
print('permutations',NP,'| arms',ARMS,flush=True)

def mcc(y,yh,thr=0.5):
    p=(yh>=thr).astype(float); tp=int(((p==1)&(y==1)).sum()); fn=int(((p==0)&(y==1)).sum())
    fp=int(((p==1)&(y==0)).sum()); tn=int(((p==0)&(y==0)).sum())
    den=np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return dict(TP=tp,FN=fn,FP=fp,TN=tn,mcc=float((tp*tn-fp*fn)/den) if den else 0.0)

for k in KEYS:
    y=Y[k]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    fns={'A1_top15' : lambda yy: ao.alt_nested_cv(D1,G,yy,1,15)[2],
         'A1_top6'  : lambda yy: ao.alt_nested_cv(D1,G,yy,1,6)[2],   # Revision: top-6 (prun_a1_k.py, seed 91) with corrected fallback branch
         'A2_top15' : lambda yy: ng.nested_cv_a2_topk(D1,yy,1,15)[2],
         'A1_viptop': lambda yy: ao.alt_viptop(D1,G,yy,1,15,2)[2],
         'A1_viptop1': lambda yy: ao.alt_viptop(D1,G,yy,1,15,1)[2],
         'A2_viptop': lambda yy: a2vip.nested_cv_a2_viptop(D1,yy,1,15,2)[2]}
    for nm in ARMS:
        p=os.path.join(CACHE,'%s_%s.json'%(k,nm))
        if os.path.exists(p):
            print('cache',k,nm,flush=True); continue
        f=fns[nm]; t=time.time()
        yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
        for j in range(NP):
            yp=np.random.default_rng([SEED[nm],j]).permutation(y); h=f(yp)
            qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
        qs=np.array(qs); as_=np.array(as_)
        res=dict(q2=oq,auc=oa,
                 p_q2=pval(qs,oq),p_auc=pval(as_,oa),
                 nperm=NP,null_q2_median=float(np.median(qs)),null_auc_mean=float(as_.mean()),
                 blowups=int(np.sum(np.abs(qs)>10)),seconds=round(time.time()-t,1),**mcc(y,yh))
        json.dump(res,open(p,'w'),indent=1,default=float)
        print('%s %-10s Q2 %6.3f p %.4f AUC %.3f p %.4f  %.0fs'%(k,nm,oq,res['p_q2'],oa,res['p_auc'],time.time()-t),flush=True)
print('done',flush=True)
