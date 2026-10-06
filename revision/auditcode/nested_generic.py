FB=[0]
def _fb():
    FB[0]+=1
    return 0.5

import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 6 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

BIG = _os.environ.get("QSAR_ROUTEA_XLSX", _os.path.join(_DATA, XLSX))
_xl=pd.ExcelFile(BIG); _clean=pd.read_excel(_xl,'Cleaned_Steps1-4')

def build(sheet, ycol=None):
    """Assemble the Route A design matrix from one sheet of Additional file 6.

    ycol names the activity column. It defaults to the third column, which is the
    convention of the strain-specific CS1-5_* sheets. The pooled sheet
    "Cleaned_Steps1-4" carries all five activity columns side by side, so ycol must
    be given explicitly there; see the strain-to-column table in the README.
    """
    hi=pd.read_excel(_xl,sheet)
    DESC=[c for c in _clean.columns if '::' in c and c in hi.columns]
    X=hi[DESC].apply(pd.to_numeric,errors='coerce').values.astype(float)
    ycol=hi.columns[2] if ycol is None else ycol; Y=(hi[ycol].astype(str).str.strip()=='active').astype(float).values
    N=len(Y); FOLDS=[]
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False; Xtr=X[tr]; R=rankdata(Xtr,axis=0)
        FOLDS.append(dict(tr=tr,i=i,Xtr=Xtr,R=R,Xte=X[i:i+1]))
    return dict(DESC=DESC,X=X,Y=Y,N=N,FOLDS=FOLDS,ycol=ycol)

def eff_filter(R,ytr,ntr):
    n1=int(ytr.sum()); n0=ntr-n1
    U1=R[ytr.astype(bool)].sum(0)-n1*(n1+1)/2.0
    rrb=2*U1/(n1*n0)-1.0
    mu=n1*n0/2.0; sigma=np.sqrt(n1*n0*(ntr+1)/12.0)
    z=(np.abs(U1-mu)-0.5)/sigma
    return np.abs(rrb), 2*norm.sf(z)

def group_reps(Xk,Rk,effk,ntr):
    k0=Xk.shape[1]
    if k0<=1: return np.arange(k0,dtype=int)
    def cor(M):
        C=np.corrcoef(M,rowvar=False); C=np.nan_to_num(C,nan=0.0); np.fill_diagonal(C,0.0)
        with np.errstate(divide='ignore',invalid='ignore'): tt=C*np.sqrt((ntr-2)/(1-C**2))
        p=2*tdist.sf(np.abs(tt),ntr-2); p[np.abs(C)>=1]=0.0
        return (np.abs(C)>0.7)&(p<0.05)
    ADJ=cor(Xk)|cor(Rk); np.fill_diagonal(ADJ,False)
    remaining=np.ones(k0,bool); reps=[]
    while remaining.any():
        idx=np.where(remaining)[0]
        deg=ADJ[np.ix_(idx,idx)].sum(1)
        if deg.max()==0: reps.extend(idx.tolist()); break
        cand=idx[deg==deg.max()]; hub=cand[np.argmax(effk[cand])]
        members=np.concatenate([[hub],np.where(ADJ[hub]&remaining)[0]])
        reps.append(int(members[np.argmax(effk[members])])); remaining[members]=False
    return np.array(sorted(reps),dtype=int)

def nested_cv(D,yv,nlv=1):
    X=D['X']; N=D['N']; yhat=np.zeros(N)
    for F in D['FOLDS']:
        tr=F['tr']; ytr=yv[tr]; ntr=int(tr.sum()); n1=int(ytr.sum()); n0=ntr-n1
        if n1<1 or n0<1: yhat[F['i']]= _fb(); continue
        eff,p=eff_filter(F['R'],ytr,ntr); keep=np.where((eff>0.3)&(p<=0.05))[0]
        if len(keep)<1: yhat[F['i']]= _fb(); continue
        rl=group_reps(F['Xtr'][:,keep],F['R'][:,keep],eff[keep],ntr); cols=keep[rl]
        Xtr=F['Xtr'][:,cols]; Xte=F['Xte'][:,cols]
        mu=Xtr.mean(0); sd=Xtr.std(0,ddof=1); sd[sd==0]=1.0
        m=PLSRegression(n_components=min(nlv,len(cols),ntr-1),scale=False).fit((Xtr-mu)/sd,ytr)
        yhat[F['i']]=m.predict((Xte-mu)/sd).ravel()[0]
    q2=1-np.sum((yv-yhat)**2)/np.sum((yv-yv.mean())**2)
    return q2,roc_auc_score(yv,yhat),yhat

def _vip_scores(model, Xs):
    import numpy as np
    T = model.x_scores_; W = model.x_weights_; Q = model.y_loadings_.ravel(); A = T.shape[1]
    ssy = np.array([(Q[a]**2)*(T[:,a]**2).sum() for a in range(A)]); p = Xs.shape[1]; V = np.zeros(p)
    for j in range(p):
        wn = (W[j,:]/np.linalg.norm(W,axis=0))**2
        V[j] = np.sqrt(p*np.sum(ssy*wn)/ssy.sum())
    return V

def _viptop(X, y, cols, k=15, nlv_vip=2):
    import numpy as np
    from sklearn.cross_decomposition import PLSRegression
    Xr = X[:, cols]
    mu = Xr.mean(0); sd = Xr.std(0, ddof=1); sd[sd == 0] = 1
    m = PLSRegression(n_components=min(nlv_vip, Xr.shape[1]), scale=False).fit((Xr-mu)/sd, y)
    V = _vip_scores(m, (Xr-mu)/sd)
    order = np.argsort(-V)[:k]
    return cols[order]

def a2_topk_select(R, ytr, ntr, keep_thr_eff=0.3, keep_thr_p=0.05, top_k=15):
    """Route A2 order (filter first, then group), further restricted to the top_k
    representatives by effect size. Returns column indices into the original 2208."""
    import numpy as np
    eff, pv = eff_filter(R, ytr, ntr)
    keep = np.where((eff > keep_thr_eff) & (pv <= keep_thr_p))[0]
    if len(keep) < 1:
        return keep, eff
    return keep, eff  # placeholder; grouping handled by caller with Xtr

def nested_cv_a2_topk(D, yv, nlv=1, top_k=15):
    import numpy as np
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.metrics import roc_auc_score
    N = D['N']; yhat = np.zeros(N)
    for F in D['FOLDS']:
        tr = F['tr']; ytr = yv[tr]; ntr = int(tr.sum())
        n1 = int(ytr.sum())
        if n1 < 1 or n1 == ntr:
            yhat[F['i']] = _fb(); continue
        eff, pv = eff_filter(F['R'], ytr, ntr)
        keep = np.where((eff > 0.3) & (pv <= 0.05))[0]
        if len(keep) < 1:
            yhat[F['i']] = _fb(); continue
        rl = group_reps(F['Xtr'][:, keep], F['R'][:, keep], eff[keep], ntr)
        reps = keep[rl]
        cols = reps[np.argsort(-eff[reps])[:top_k]] if len(reps) > top_k else reps
        Xtr = F['Xtr'][:, cols]; Xte = F['Xte'][:, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(nlv, len(cols), ntr-1), scale=False).fit((Xtr-mu)/sd, ytr)
        yhat[F['i']] = m.predict((Xte-mu)/sd).ravel()[0]
    q2v = 1 - np.sum((yv-yhat)**2) / np.sum((yv-yv.mean())**2)
    return q2v, roc_auc_score(yv, yhat), yhat
