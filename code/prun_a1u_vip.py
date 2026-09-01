"""VIP selection starting from the label-free Route A1u descriptor set.
A1u selection -> first PLS-DA (1 LV) -> top-15 by VIP -> second PLS-DA (1 LV)."""
import sys, os, time, numpy as np, nested_generic as ng, altorder as ao, unsup_rep as ur
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

def viptop_unsup(D, reps_u, yv, nlv=1, top_k=15, nlv_vip=1):
    N=D['N']; yhat=np.zeros(N)
    for F in D['FOLDS']:
        tr=F['tr']; ytr=yv[tr]; ntr=int(tr.sum()); n1=int(ytr.sum())
        if n1<1 or n1==ntr: yhat[F['i']]=ytr.mean(); continue
        cols,eff,pv = ur.unsup_select(reps_u, F['R'], ytr, ntr)
        if len(cols)<1: yhat[F['i']]=ytr.mean(); continue
        if len(cols)>top_k: cols=ng._viptop(F['Xtr'], ytr, cols, top_k, nlv_vip)
        Xtr=F['Xtr'][:,cols]; Xte=F['Xte'][:,cols]
        mu=Xtr.mean(0); sd=Xtr.std(0,ddof=1); sd[sd==0]=1.0
        m=PLSRegression(n_components=min(nlv,len(cols),ntr-1),scale=False).fit((Xtr-mu)/sd,ytr)
        yhat[F['i']]=m.predict((Xte-mu)/sd).ravel()[0]
    q2=1-np.sum((yv-yhat)**2)/np.sum((yv-yv.mean())**2)
    return q2, roc_auc_score(yv,yhat), yhat

if __name__=='__main__':
    short=sys.argv[1]; sheet=sys.argv[2]; target=int(sys.argv[3]); code=sys.argv[4]
    budget=float(sys.argv[5]) if len(sys.argv)>5 else 250.0
    STORE=f'/home/claude/perm_a1uvip_{code}.npz'; t0=time.time()
    D=ng.build(sheet); groups=ao.get_groups(D['X']); reps=ur.centrality_reps(D['X'],groups)
    if os.path.exists(STORE):
        d=np.load(STORE); q2s=list(d['q2s']); aucs=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
    else:
        oq,oa,_=viptop_unsup(D,reps,D['Y']); q2s=[]; aucs=[]
    def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
    while len(q2s)<target and (time.time()-t0)<budget:
        k=len(q2s); yv=np.random.default_rng([313,k]).permutation(D['Y'])
        q,a,_=viptop_unsup(D,reps,yv); q2s.append(q); aucs.append(a)
        if len(q2s)%200==0: save()
    save(); n=len(q2s)
    pQ=(1+np.sum(np.array(q2s)>=oq))/(1+n); pA=(1+np.sum(np.array(aucs)>=oa))/(1+n)
    print(f'{short} A1u-VIP15: n={n} | Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | {time.time()-t0:.0f}s')
