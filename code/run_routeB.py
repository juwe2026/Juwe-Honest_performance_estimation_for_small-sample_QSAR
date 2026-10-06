# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
"""Route B permutation test with a FIXED, a priori descriptor set (leakage-free:
apparent == nested, and a simple permutation test is already honest)."""
import sys, os, pickle, time, numpy as np, pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

short=sys.argv[1]; which=sys.argv[2]; target=int(sys.argv[3])
budget=float(sys.argv[4]) if len(sys.argv)>4 else 230.0
STORE=f'/home/claude/permB_{which}_{short.replace(".","")}.npz'
t0=time.time()
full=pd.read_pickle('/home/claude/desc_full.pkl')
cols = list(np.load('/home/claude/lit25_desc.npy',allow_pickle=True)) if which=='25' \
       else list(np.load('/home/claude/reduced_desc.npy',allow_pickle=True))
Yc={'H.I.':'H. I. (H. influenzae)','S.A.':'S. A. (Staph. Aureus)','S.Pneu':'S. Pneu (S. pneumoniae)',
    'S.Pyo':'S. pyo (S. pyogenes)','P.A.':'P. A. (P. aeruginosa)'}
y=(full[Yc[short]].astype(str).str.strip().str.lower()=='active').astype(float).values
X=full[cols].astype(float).values
N=len(y)

def loo(X,y,nlv=1):
    yh=np.zeros(N)
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False
        Xt=X[tr]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1); sd[sd==0]=1
        m=PLSRegression(n_components=min(nlv,X.shape[1],N-2),scale=False).fit((Xt-mu)/sd,y[tr])
        yh[i]=m.predict((X[i:i+1]-mu)/sd).ravel()[0]
    return yh
def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)

if os.path.exists(STORE):
    d=np.load(STORE); qs=list(d['q2s']); as_=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
else:
    yh=loo(X,y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[]; as_=[]
def save(): np.savez(STORE,q2s=np.array(qs),aucs=np.array(as_),obs_q2=oq,obs_auc=oa)
while len(qs)<target and (time.time()-t0)<budget:
    k=len(qs); yp=np.random.default_rng([701 if which=='25' else 702,k]).permutation(y)
    yhp=loo(X,yp); qs.append(q2(yp,yhp)); as_.append(roc_auc_score(yp,yhp))
    if len(qs)%200==0: save()
save(); n=len(qs)
pq=(1+np.sum(np.array(qs)>=oq - 1e-9))/(1+n); pa=(1+np.sum(np.array(as_)>=oa - 1e-9))/(1+n)
print(f'{short} RouteB-{which}: n={n} | Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pq:.4f} p(AUC)={pa:.4f} | {time.time()-t0:.0f}s')
