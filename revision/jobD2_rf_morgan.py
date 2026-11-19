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
import os, json
# G1-3: random forest control (min_samples_leaf=1), OOB estimate observed + 500 label permutations (seed [8501,j])
from common import *
import core
from sklearn.ensemble import RandomForestClassifier
MB=np.load((_OUT + '/morgan2_bits.npy'))
SETS={'Morgan2_bits':MB}
NP=500
def rf_oob(Xm,y,seed=0):
    rf=RandomForestClassifier(300,min_samples_leaf=1,max_features='sqrt',oob_score=True,n_jobs=1,random_state=seed).fit(Xm,y)
    return rf.oob_decision_function_[:,1]
def rf_loo(Xm,y):
    yh=np.zeros(N)
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False
        rf=RandomForestClassifier(300,min_samples_leaf=1,max_features='sqrt',n_jobs=1,random_state=0).fit(Xm[tr],y[tr]); yh[i]=rf.predict_proba(Xm[i:i+1])[0,1]
    return yh
out={}
_pp=(_OUT + '/D_rf_morgan.json')
if os.environ.get('RESUME')=='1' and os.path.exists(_pp): out=json.load(open(_pp))   # restart after a change of container
for nm,Xm in SETS.items():
    out.setdefault(nm,{})
    for k in KEYS:
        if k in out[nm]: print('cache',nm,k,flush=True); continue
        t=time.time(); y=Y[k]; po=rf_oob(Xm,y); oq=q2(y,po); oa=roc_auc_score(y,po)
        pl=rf_loo(Xm,y)
        qs=[];as_=[]
        for j in range(NP):
            yp=np.random.default_rng([8501,j]).permutation(y); p=rf_oob(Xm,yp); qs.append(q2(yp,p)); as_.append(roc_auc_score(yp,p))
        out[nm][k]=dict(oob_q2=oq,oob_auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),loo_q2=q2(y,pl),loo_auc=roc_auc_score(y,pl),loo=mcc_conf(y,pl),nperm=NP)
        print(nm,k,'OOB Q2 %.3f p %.3f AUC %.3f p %.3f | LOO Q2 %.3f AUC %.3f MCC %.3f | %.0fs'%(oq,out[nm][k]['p_q2'],oa,out[nm][k]['p_auc'],q2(y,pl),roc_auc_score(y,pl),out[nm][k]['loo']['mcc'],time.time()-t),flush=True)
        dump('D_rf_morgan',out)
