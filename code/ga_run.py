"""GA-Deskriptorauswahl nach Mansouri-Bauplan, apparent und vollstaendig genestet."""
import sys, os, time, pickle, numpy as np
import nested_generic as ng, altorder as ao, unsup_rep as ur
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

def q2(a,b): return 1-np.sum((a-b)**2)/np.sum((a-a.mean())**2)

def fitness(chrom, Xs, ys, pen=0.002):
    idx=np.where(chrom)[0]
    if len(idx)<2: return -9.9
    k=int(min(5, ys.sum(), (ys==0).sum()))
    if k<2: return -9.9
    skf=StratifiedKFold(n_splits=k, shuffle=True, random_state=0); yh=np.zeros(len(ys))
    for tr,te in skf.split(Xs,ys):
        A=Xs[tr][:,idx]; mu=A.mean(0); sd=A.std(0,ddof=1); sd[sd==0]=1
        m=PLSRegression(1,scale=False).fit((A-mu)/sd,ys[tr])
        yh[te]=m.predict((Xs[te][:,idx]-mu)/sd).ravel()
    return q2(ys,yh)-pen*len(idx)

def ga(Xs, ys, npop=30, ngen=50, pc=0.5, pm=0.01, seed=0):
    rng=np.random.default_rng(seed); P=Xs.shape[1]
    pop=rng.random((npop,P))<0.05; fit=np.array([fitness(c,Xs,ys) for c in pop])
    for g in range(ngen):
        o=np.argsort(-fit); pop=pop[o]; fit=fit[o]; new=pop[:npop//2].copy()
        while len(new)<npop:
            a,b=rng.integers(0,npop//2,2); cut=rng.integers(1,P)
            ch=np.concatenate([pop[a][:cut],pop[b][cut:]]) if rng.random()<pc else pop[a].copy()
            ch=ch^(rng.random(P)<pm); new=np.vstack([new,ch])
        pop=new; fit=np.array([fitness(c,Xs,ys) for c in pop])
    b=int(np.argmax(fit)); return pop[b], fit[b]

def loo_fixed(Xs, ys, idx):
    n=len(ys); yh=np.zeros(n)
    for i in range(n):
        tr=np.ones(n,bool); tr[i]=False
        A=Xs[tr][:,idx]; mu=A.mean(0); sd=A.std(0,ddof=1); sd[sd==0]=1
        m=PLSRegression(1,scale=False).fit((A-mu)/sd,ys[tr])
        yh[i]=m.predict((Xs[i:i+1,idx]-mu)/sd).ravel()[0]
    return yh

if __name__=='__main__':
    code=sys.argv[1]; sheet=sys.argv[2]; budget=float(sys.argv[3]) if len(sys.argv)>3 else 250.0
    STORE=f'/home/claude/ga_{code}.pkl'; t0=time.time()
    D=ng.build(sheet); y=D['Y']; N=D['N']
    D0=ng.build('CS1-4_H.I.'); groups=ao.get_groups(D0['X']); reps=ur.centrality_reps(D0['X'],groups)
    X=D['X'][:,reps]
    S=pickle.load(open(STORE,'rb')) if os.path.exists(STORE) else {'yhn':np.zeros(N),'done':[],'sizes':[]}
    if 'app' not in S:
        ch,f=ga(X,y); idx=np.where(ch)[0]
        yh=loo_fixed(X,y,idx)
        S['app']=dict(n=len(idx), fit_q2=float(f+0.002*len(idx)), loo_q2=float(q2(y,yh)), loo_auc=float(roc_auc_score(y,yh)))
        pickle.dump(S,open(STORE,'wb'))
    for i in range(N):
        if i in S['done']: continue
        if time.time()-t0>budget: break
        tr=np.ones(N,bool); tr[i]=False
        c2,_=ga(X[tr],y[tr],seed=100+i); j2=np.where(c2)[0]
        A=X[tr][:,j2]; mu=A.mean(0); sd=A.std(0,ddof=1); sd[sd==0]=1
        m=PLSRegression(1,scale=False).fit((A-mu)/sd,y[tr])
        S['yhn'][i]=m.predict((X[i:i+1,j2]-mu)/sd).ravel()[0]
        S['done'].append(i); S['sizes'].append(int(len(j2))); pickle.dump(S,open(STORE,'wb'))
    d=len(S['done'])
    print(f"{code}: {d}/{N} Folds fertig ({time.time()-t0:.0f}s)")
    if d==N:
        print(f"  apparent: GA-Fitness Q2={S['app']['fit_q2']:.3f} ({S['app']['n']} Desk.) | LOO fix Q2={S['app']['loo_q2']:.3f} AUC={S['app']['loo_auc']:.3f}")
        print(f"  genestet: Q2={q2(y,S['yhn']):.3f} AUC={roc_auc_score(y,S['yhn']):.3f} | Modellgroesse Median {int(np.median(S['sizes']))}, {min(S['sizes'])}-{max(S['sizes'])}")
