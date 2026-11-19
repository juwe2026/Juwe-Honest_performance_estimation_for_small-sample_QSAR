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
import sys,os,pickle,json,time
# The revision scripts were run on the Route B descriptor file as submitted. QSAR_ROUTEB_XLSX overrides it.
os.environ.setdefault('QSAR_ROUTEB_XLSX', _ROOT + '/data/RouteB_descriptor_data_submitted.xlsx')
sys.path.insert(0,(_ROOT + '/code'))
os.environ['QSAR_RECALC_DIR']=_WORK
import numpy as np
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
B=pickle.load(open((_WORK + '/base.pkl'),'rb'))
D=B['D']; G=B['g']; REPS_U=np.asarray(B['reps']); ADJ=B['ADJ']
X=D['HI']['X']; DESC=D['HI']['DESC']; N=31
KEYS=['HI','SA','SPneu','SPyo','PA']
Y={k:D[k]['Y'] for k in KEYS}
def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)
def loo_pls(Xm,y,nlv=1):
    yh=np.zeros(len(y))
    for i in range(len(y)):
        tr=np.ones(len(y),bool); tr[i]=False
        Xt=Xm[tr]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1)
        ok=sd>1e-6   # drop columns that are numerically constant (SD < 1e-6) in the training fold
        Xt=Xt[:,ok]; mu=mu[ok]; sd=sd[ok]
        m=PLSRegression(n_components=min(nlv,Xt.shape[1],len(y)-2),scale=False).fit((Xt-mu)/sd,y[tr])
        yh[i]=m.predict((Xm[i:i+1][:,ok]-mu)/sd).ravel()[0]
    return yh
exec(open((_REV + '/pval_tie.py')).read())   # treatment of ties, revision of 03.10.2026
def mcc_conf(y,yh,thr=0.5):
    p=(yh>=thr).astype(float); tp=int(((p==1)&(y==1)).sum()); fn=int(((p==0)&(y==1)).sum())
    fp=int(((p==1)&(y==0)).sum()); tn=int(((p==0)&(y==0)).sum())
    den=np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return dict(TP=tp,FN=fn,FP=fp,TN=tn,sens=tp/(tp+fn),spec=tn/(tn+fp),mcc=float((tp*tn-fp*fn)/den) if den else 0.0)
def dump(name,obj): json.dump(obj,open((_OUT + '/%s.json')%name,'w'),indent=1,default=float)
