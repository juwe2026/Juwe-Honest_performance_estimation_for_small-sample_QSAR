import numpy as np, pickle, nested_generic as ng, altorder as ao
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression

sheets={'H.I.':'CS1-4_H.I.','S.A.':'CS1-4_S.A.','S.Pneu':'CS1-4_S. Pneu','S.Pyo':'CS1-4_S. Pyo','P.A.':'CS1-4_P.A.'}
order=['H.I.','S.A.','S.Pneu','S.Pyo','P.A.']

def vip_scores(model,Xs):
    T=model.x_scores_; W=model.x_weights_; Q=model.y_loadings_.ravel(); A=T.shape[1]
    ssy=np.array([(Q[a]**2)*(T[:,a]**2).sum() for a in range(A)]); p=Xs.shape[1]; V=np.zeros(p)
    for j in range(p):
        wn=(W[j,:]/np.linalg.norm(W,axis=0))**2; V[j]=np.sqrt(p*np.sum(ssy*wn)/ssy.sum())
    return V

def fit_vip(X,y,nlv=2):
    mu=X.mean(0); sd=X.std(0,ddof=1); sd[sd==0]=1; Xs=(X-mu)/sd
    m=PLSRegression(min(nlv,Xs.shape[1]),scale=False).fit(Xs,y)
    return vip_scores(m,Xs), m

D0=ng.build('CS1-4_H.I.'); groups=ao.get_groups(D0['X'])
DATA={}
for s,sh in sheets.items():
    D=ng.build(sh); y=D['Y']; N=D['N']; R=rankdata(D['X'],axis=0)
    # --- A1 top-15 by effect size (apparent, on all data) ---
    cols,eff,pv,reps=ao.alt_select(groups,R,y,N)
    cols15=cols[np.argsort(-eff[cols])[:15]] if len(cols)>15 else cols
    names15=[D['DESC'][c] for c in cols15]
    Va,mfit=fit_vip(D['X'][:,cols15], y, 2)          # apparent VIP
    effA={D['DESC'][c]:eff[c] for c in cols15}

    # --- mean nested VIP: repeat A1 selection (top-15 by eff.size) inside each LOO fold ---
    accV={n:[] for n in names15}; freq={n:0 for n in names15}
    for F in D['FOLDS']:
        tr=F['tr']; ytr=y[tr]; ntr=int(tr.sum())
        cF,effF,pvF,_=ao.alt_select(groups,F['R'],ytr,ntr)
        if len(cF)<1: continue
        cF15=cF[np.argsort(-effF[cF])[:15]] if len(cF)>15 else cF
        Vf,_=fit_vip(F['Xtr'][:,cF15], ytr, 2)
        for j,c in enumerate(cF15):
            nm=D['DESC'][c]
            if nm in accV: accV[nm].append(Vf[j]); freq[nm]+=1

    rowsA=[]
    for j,c in enumerate(cols15):
        nm=D['DESC'][c]
        rowsA.append(dict(Descriptor=nm, VIP_apparent=round(float(Va[j]),3),
                          VIP_nested_mean=round(float(np.mean(accV[nm])),3) if accV[nm] else None,
                          Sel_freq=f"{freq[nm]}/31", EffectSize=round(float(effA[nm]),3)))
    DATA[s]=dict(chemdes=rowsA, names15=names15,
                 scoreA=dict(X=D['X'][:,cols15].copy(), y=y.copy(), names=names15))
    print(f'{s}: A1-top15 recomputed ({len(rowsA)} descriptors)')

pickle.dump(DATA, open('/home/claude/vip_data_A1.pkl','wb'))
print('gespeichert: vip_data_A1.pkl')
