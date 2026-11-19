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
import numpy as np, pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
full=pd.read_pickle(_WORK + '/desc_full.pkl'); reduced=list(np.load(_WORK + '/reduced_desc.npy',allow_pickle=True))
ACT=[c for c in full.columns if c.startswith(('H. I','P. A','S. A','S. Pneu','S. pyo'))]
NAME={'H. I. (H. influenzae)':'H.I.','P. A. (P. aeruginosa)':'P.A.','S. A. (Staph. Aureus)':'S.A.','S. Pneu (S. pneumoniae)':'S.Pneu','S. pyo (S. pyogenes)':'S.Pyo'}
X=full[reduced].astype(float).values
def loo(X,y,nlv):
    n=len(y); yh=np.zeros(n)
    for i in range(n):
        tr=np.ones(n,bool);tr[i]=False
        Xtr,Xte=X[tr],X[i:i+1]; mu=Xtr.mean(0);sd=Xtr.std(0,ddof=1);sd[sd==0]=1
        m=PLSRegression(min(nlv,tr.sum()-1),scale=False).fit((Xtr-mu)/sd,y[tr]); yh[i]=m.predict((Xte-mu)/sd).ravel()[0]
    return yh
def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)
rng=np.random.default_rng(7)
print(f"{'strain':7} {'nLV=1: Q2':>10} {'AUC':>6} {'p(Q2)':>7} {'p(AUC)':>7}")
res={}
for col in ACT:
    y=(full[col].astype(str).str.strip().str.lower()=='active').astype(float).values; short=NAME[col]
    yh=loo(X,y,1); Q=q2(y,yh); A=roc_auc_score(y,yh)
    pq=pa=np.empty(2000)
    PQ=np.empty(2000);PA=np.empty(2000)
    for k in range(2000):
        yp=rng.permutation(y); yhp=loo(X,yp,1); PQ[k]=q2(yp,yhp); PA[k]=roc_auc_score(yp,yhp)
    pQ=(1+np.sum(PQ>=Q - 1e-9))/2001; pA=(1+np.sum(PA>=A - 1e-9))/2001
    res[short]=(Q,A,pQ,pA)
    print(f'{short:7} {Q:10.3f} {A:6.3f} {pQ:7.4f} {pA:7.4f}')
np.save((_WORK + '/lit_nlv1.npy'),res)
