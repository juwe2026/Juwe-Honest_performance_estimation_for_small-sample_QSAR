import numpy as np, pandas as pd, openpyxl
from openpyxl.styles import Font, PatternFill
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score, matthews_corrcoef

full=pd.read_pickle('desc_full.pkl')
reduced=list(np.load('reduced_desc.npy',allow_pickle=True))
ACT=[c for c in full.columns if c.startswith(('H. I','P. A','S. A','S. Pneu','S. pyo'))]
ALLD=[c for c in full.columns if c not in ['Monoterpenoid','SMILES']+ACT]
NAME={'H. I. (H. influenzae)':'H.I.','P. A. (P. aeruginosa)':'P.A.','S. A. (Staph. Aureus)':'S.A.',
      'S. Pneu (S. pneumoniae)':'S.Pneu','S. pyo (S. pyogenes)':'S.Pyo'}

def asc(Xtr,Xte):
    mu=Xtr.mean(0); sd=Xtr.std(0,ddof=1); sd[sd==0]=1; return (Xtr-mu)/sd,(Xte-mu)/sd
def loo(X,y,nlv,selector=None):
    n=len(y); yh=np.zeros(n)
    for i in range(n):
        tr=np.ones(n,bool); tr[i]=False
        Xtr,Xte=X[tr],X[i:i+1]; ytr=y[tr]
        cols=selector(Xtr,ytr,nlv) if selector else np.arange(X.shape[1])
        if len(cols)<1: yh[i]=ytr.mean(); continue
        a,b=asc(Xtr[:,cols],Xte[:,cols]); use=min(nlv,len(cols),tr.sum()-1)
        m=PLSRegression(use,scale=False).fit(a,ytr); yh[i]=m.predict(b).ravel()[0]
    return yh
def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)
def vip(model,X,y):
    T=model.x_scores_; W=model.x_weights_; Q=model.y_loadings_.ravel(); A=T.shape[1]
    ssy=np.array([(Q[a]**2)*(T[:,a]**2).sum() for a in range(A)]); p=X.shape[1]; V=np.zeros(p)
    for j in range(p):
        wn=(W[j,:]/np.linalg.norm(W,axis=0))**2; V[j]=np.sqrt(p*np.sum(ssy*wn)/ssy.sum())
    return V
def fit_vip(X,y,nlv):
    mu=X.mean(0); sd=X.std(0,ddof=1); sd[sd==0]=1; m=PLSRegression(nlv,scale=False).fit((X-mu)/sd,y)
    return vip(m,(X-mu)/sd,y)
def metrics(y,yh):
    p=(yh>=0.5).astype(int); TP=int(((p==1)&(y==1)).sum());FN=int(((p==0)&(y==1)).sum())
    FP=int(((p==1)&(y==0)).sum());TN=int(((p==0)&(y==0)).sum())
    sens=TP/(TP+FN) if TP+FN else np.nan; spec=TN/(TN+FP) if TN+FP else np.nan
    try: mcc=matthews_corrcoef(y,p)
    except: mcc=np.nan
    return dict(TP=TP,FN=FN,FP=FP,TN=TN,sens=sens,spec=spec,ba=np.nanmean([sens,spec]),
                acc=(TP+TN)/len(y),mcc=mcc)

Xred=full[reduced].astype(float).values
Xall=full[ALLD].astype(float).values
rng=np.random.default_rng(2026)

rows=[]; vip_tab={'Deskriptor':reduced}; twostage=[]
for col in ACT:
    y=(full[col].astype(str).str.strip().str.lower()=='active').astype(float).values
    short=NAME[col]; n1=int(y.sum())
    # --- Hauptmodell: reduzierter Satz (unsupervised reduziert -> einfache Permutation ist ehrlich) ---
    scan=[(k,q2(y,loo(Xred,y,k))) for k in range(1,5)]
    nlv=sorted(scan,key=lambda r:(-r[1],r[0]))[0][0]
    yh=loo(Xred,y,nlv); Q=q2(y,yh); AUC=roc_auc_score(y,yh); mt=metrics(y,yh)
    # ehrliche Permutation (fixer Satz)
    pq=pa=0
    permsQ=np.empty(2000); permsA=np.empty(2000)
    for k in range(2000):
        yp=rng.permutation(y); yhp=loo(Xred,yp,nlv)
        permsQ[k]=q2(yp,yhp); permsA[k]=roc_auc_score(yp,yhp)
    pQ=(1+np.sum(permsQ>=Q))/2001; pA=(1+np.sum(permsA>=AUC))/2001
    # VIP auf Gesamtdaten (nur Interpretation)
    mu=Xred.mean(0); sd=Xred.std(0,ddof=1); sd[sd==0]=1
    Vfull=vip(PLSRegression(nlv,scale=False).fit((Xred-mu)/sd,y),(Xred-mu)/sd,y)
    vip_tab[short]=np.round(Vfull,3)
    rows.append(dict(Erreger=short,aktiv=n1,nLV=nlv,R2Y=None,Q2=round(Q,3),AUC=round(AUC,3),
                     BA=round(mt['ba'],3),Sens=round(mt['sens'],3),Spez=round(mt['spec'],3),
                     MCC=round(mt['mcc'],3),Konfusion=f"{mt['TP']}/{mt['FN']}/{mt['FP']}/{mt['TN']}",
                     p_Q2=round(pQ,4),p_AUC=round(pA,4)))
    # --- Zweistufiger VIP-Ansatz (Start: volle 25) : naiv vs. genestet ---
    def sel_vip_train(Xtr,ytr,nlv):  # Selektion NUR auf Training (genestet)
        V=fit_vip(Xtr,ytr,min(nlv,Xtr.shape[1])); c=np.where(V>1)[0]; return c if len(c) else np.arange(Xtr.shape[1])
    # naiv: VIP>1 auf ALLEN Daten bestimmen, dann fixe Spalten in LOO
    Vall=fit_vip(Xall,y,nlv); naive_cols=np.where(Vall>1)[0]
    yh_naive=loo(Xall[:,naive_cols] if len(naive_cols) else Xall, y, nlv)
    Qn=q2(y,yh_naive); An=roc_auc_score(y,yh_naive)
    # genestet: VIP>1 in jeder Fold neu
    yh_nest=loo(Xall,y,nlv,selector=sel_vip_train); Qg=q2(y,yh_nest); Ag=roc_auc_score(y,yh_nest)
    twostage.append(dict(Erreger=short,aktiv=n1,nLV=nlv,
                         Q2_naiv=round(Qn,3),AUC_naiv=round(An,3),
                         Q2_genestet=round(Qg,3),AUC_genestet=round(Ag,3),
                         nVIP_gt1_gesamt=int(len(naive_cols))))
    print(f'{short}: nLV={nlv} | reduziert Q2={Q:.3f} AUC={AUC:.3f} p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | VIP naiv Q2={Qn:.3f} genestet Q2={Qg:.3f}')

RES=pd.DataFrame(rows); VIP=pd.DataFrame(vip_tab); TWO=pd.DataFrame(twostage)
fn='/mnt/user-data/outputs/20260706_31MT_LitDeskriptoren_Analyse.xlsx'
with pd.ExcelWriter(fn,engine='openpyxl',mode='a',if_sheet_exists='replace') as w:
    RES.to_excel(w,sheet_name='PLS-DA_reduziert',index=False)
    VIP.to_excel(w,sheet_name='VIP_reduziert',index=False)
    TWO.to_excel(w,sheet_name='VIP_zweistufig',index=False)
RES.to_pickle('/home/claude/lit_results.pkl'); TWO.to_pickle('/home/claude/lit_twostage.pkl')
print('\nGeschrieben.')
