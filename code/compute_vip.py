import numpy as np, pandas as pd, pickle
from scipy.stats import rankdata, norm
from sklearn.cross_decomposition import PLSRegression
import nested_generic as ng

def vip(model,X):
    T=model.x_scores_; W=model.x_weights_; Q=model.y_loadings_.ravel(); A=T.shape[1]
    ssy=np.array([(Q[a]**2)*(T[:,a]**2).sum() for a in range(A)]); p=X.shape[1]; V=np.zeros(p)
    for j in range(p):
        wn=(W[j,:]/np.linalg.norm(W,axis=0))**2; V[j]=np.sqrt(p*np.sum(ssy*wn)/ssy.sum())
    return V
def fit_vip(X,y,nlv=2):
    mu=X.mean(0);sd=X.std(0,ddof=1);sd[sd==0]=1; Xs=(X-mu)/sd
    m=PLSRegression(min(nlv,Xs.shape[1]),scale=False).fit(Xs,y); return vip(m,Xs), (mu,sd)
def eff_rb(x,y):  # rank-biserial effect size, active vs non-active
    a=x[y==1]; b=x[y==0]
    if len(np.unique(np.concatenate([a,b])))<2: return 0.0
    r=rankdata(np.concatenate([a,b])); R1=r[:len(a)].sum(); U1=R1-len(a)*(len(a)+1)/2
    return abs(2*U1/(len(a)*len(b))-1)

sheets={'H.I.':'CS1-4_H.I.','S.A.':'CS1-4_S.A.','S.Pneu':'CS1-4_S. Pneu','S.Pyo':'CS1-4_S. Pyo','P.A.':'CS1-4_P.A.'}
# Literature data
full=pd.read_pickle('desc_full.pkl'); reduced=list(np.load('reduced_desc.npy',allow_pickle=True))
LITALL=[c for c in full.columns if c not in ['Monoterpenoid','SMILES'] and not c.startswith(('H. I','P. A','S. A','S. Pneu','S. pyo'))]
Ycols={'H.I.':'H. I. (H. influenzae)','S.A.':'S. A. (Staph. Aureus)','S.Pneu':'S. Pneu (S. pneumoniae)','S.Pyo':'S. pyo (S. pyogenes)','P.A.':'P. A. (P. aeruginosa)'}

DATA={}
for short,sh in sheets.items():
    D=ng.build(sh); y=D['Y']; N=D['N']
    # ---- Route A top-15 (ChemDes) ----
    cols_all,eff_all=ng.full_select(D); cols15=cols_all[:15]
    names15=[D['DESC'][c] for c in cols15]
    Va,_=fit_vip(D['X'][:,cols15], y, 2)                 # apparent VIP
    effA={D['DESC'][c]:eff_all[c] for c in cols15}
    # nested: per fold select top15, VIP
    accV={n:[] for n in D['DESC']}; freq={n:0 for n in D['DESC']}
    for F in D['FOLDS']:
        tr=F['tr']; ytr=y[tr]; ntr=int(tr.sum())
        e,p=ng.eff_filter(F['R'],ytr,ntr); keep=np.where((e>0.3)&(p<=0.05))[0]
        if len(keep)<1: continue
        rl=ng.group_reps(F['Xtr'][:,keep],F['R'][:,keep],e[keep],ntr); cc=keep[rl]
        cc=cc[np.argsort(-e[cc])[:15]]
        Vf,_=fit_vip(F['Xtr'][:,cc], ytr, 2)
        for j,c in enumerate(cc): nm=D['DESC'][c]; accV[nm].append(Vf[j]); freq[nm]+=1
    rowsA=[]
    for j,c in enumerate(cols15):
        nm=D['DESC'][c]
        rowsA.append(dict(Descriptor=nm, VIP_apparent=round(float(Va[j]),3),
                          VIP_nested_mean=round(float(np.mean(accV[nm])),3) if accV[nm] else None,
                          Sel_freq=f"{freq[nm]}/31", EffectSize=round(float(effA[nm]),3)))
    # ---- Route B (13 literature) ----
    ycol=Ycols[short]; yB=(full[ycol].astype(str).str.strip().str.lower()=='active').astype(float).values
    XB=full[reduced].astype(float).values
    Vb,_=fit_vip(XB,yB,2)
    accB={d:[] for d in reduced}
    for i in range(N):
        trm=np.ones(N,bool); trm[i]=False
        Vf,_=fit_vip(XB[trm],yB[trm],2)
        for j,d in enumerate(reduced): accB[d].append(Vf[j])
    effLit={d:eff_rb(full[d].astype(float).values,yB) for d in LITALL}
    # ---- Two-stage VIP (25 literature) ----
    X25=full[LITALL].astype(float).values
    V25,_=fit_vip(X25,yB,2)
    acc25={d:[] for d in LITALL}; freq25={d:0 for d in LITALL}
    for i in range(N):
        trm=np.ones(N,bool); trm[i]=False
        Vf,_=fit_vip(X25[trm],yB[trm],2)
        for j,d in enumerate(LITALL):
            acc25[d].append(Vf[j])
            if Vf[j]>1: freq25[d]+=1
    rowsLit=[]
    for j,d in enumerate(LITALL):
        inB = d in reduced
        rowsLit.append(dict(Descriptor=d,
            RouteB_VIP_apparent=round(float(Vb[reduced.index(d)]),3) if inB else None,
            RouteB_VIP_nested_mean=round(float(np.mean(accB[d])),3) if inB else None,
            TwoStage_VIP_apparent=round(float(V25[j]),3),
            TwoStage_VIP_nested_mean=round(float(np.mean(acc25[d])),3),
            TwoStage_selfreq_VIPgt1=f"{freq25[d]}/31",
            EffectSize=round(float(effLit[d]),3)))
    DATA[short]=dict(chemdes=rowsA, lit=rowsLit, names15=names15, reduced=reduced,
                     # for score plots (Route A top-15, nLV=2)
                     scoreA=dict(X=D['X'][:,cols15].copy(), y=y.copy(), names=names15))
    print(f'{short}: Route A top15 VIP ok ({len(rowsA)}), Route B/Two-stage ok ({len(rowsLit)})')

pickle.dump(DATA, open('/home/claude/vip_data.pkl','wb'))
print('gespeichert.')
