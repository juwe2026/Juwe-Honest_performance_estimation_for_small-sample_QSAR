# Revision 03.10.2026: permutation p-values count ties (tolerance 1e-9), see CHANGELOG_revision.md
# REVISION OCTOBER 2026: the empty-selection / single-class fallback below previously
# predicted ytr.mean(). Under leave-one-out that equals (S - y_i)/(n - 1) and therefore
# carries the held-out label; it is now the neutral constant 0.5. See CHANGELOG_revision.md.
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
import numpy as np, pandas as pd, openpyxl
from openpyxl.styles import Font, PatternFill
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score, matthews_corrcoef

full=pd.read_pickle(_WORK + '/desc_full.pkl')
reduced=list(np.load(_WORK + '/reduced_desc.npy',allow_pickle=True))
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
        if len(cols)<1: yh[i]=0.5; continue
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
    # --- main model: reduced set (reduced unsupervised -> simple permutation is honest) ---
    scan=[(k,q2(y,loo(Xred,y,k))) for k in range(1,5)]
    nlv=sorted(scan,key=lambda r:(-r[1],r[0]))[0][0]
    yh=loo(Xred,y,nlv); Q=q2(y,yh); AUC=roc_auc_score(y,yh); mt=metrics(y,yh)
    # honest permutation (fixed set)
    pq=pa=0
    permsQ=np.empty(2000); permsA=np.empty(2000)
    for k in range(2000):
        yp=rng.permutation(y); yhp=loo(Xred,yp,nlv)
        permsQ[k]=q2(yp,yhp); permsA[k]=roc_auc_score(yp,yhp)
    pQ=(1+np.sum(permsQ>=Q - 1e-9))/2001; pA=(1+np.sum(permsA>=AUC - 1e-9))/2001
    # VIP on all data (interpretation only)
    mu=Xred.mean(0); sd=Xred.std(0,ddof=1); sd[sd==0]=1
    Vfull=vip(PLSRegression(nlv,scale=False).fit((Xred-mu)/sd,y),(Xred-mu)/sd,y)
    vip_tab[short]=np.round(Vfull,3)
    rows.append(dict(Erreger=short,aktiv=n1,nLV=nlv,R2Y=None,Q2=round(Q,3),AUC=round(AUC,3),
                     BA=round(mt['ba'],3),Sens=round(mt['sens'],3),Spez=round(mt['spec'],3),
                     MCC=round(mt['mcc'],3),Konfusion=f"{mt['TP']}/{mt['FN']}/{mt['FP']}/{mt['TN']}",
                     p_Q2=round(pQ,4),p_AUC=round(pA,4)))
    # --- two-stage VIP approach (start: all 25): naive vs. nested ---
    def sel_vip_train(Xtr,ytr,nlv):  # selection ONLY on training (nested)
        V=fit_vip(Xtr,ytr,min(nlv,Xtr.shape[1])); c=np.where(V>1)[0]; return c if len(c) else np.arange(Xtr.shape[1])
    # naive: determine VIP>1 on ALL data, then fixed columns in LOO
    Vall=fit_vip(Xall,y,nlv); naive_cols=np.where(Vall>1)[0]
    yh_naive=loo(Xall[:,naive_cols] if len(naive_cols) else Xall, y, nlv)
    Qn=q2(y,yh_naive); An=roc_auc_score(y,yh_naive)
    # nested: VIP>1 again in every fold
    yh_nest=loo(Xall,y,nlv,selector=sel_vip_train); Qg=q2(y,yh_nest); Ag=roc_auc_score(y,yh_nest)
    twostage.append(dict(Erreger=short,aktiv=n1,nLV=nlv,
                         Q2_naiv=round(Qn,3),AUC_naiv=round(An,3),
                         Q2_genestet=round(Qg,3),AUC_genestet=round(Ag,3),
                         nVIP_gt1_gesamt=int(len(naive_cols))))
    print(f'{short}: nLV={nlv} | reduced Q2={Q:.3f} AUC={AUC:.3f} p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | VIP naive Q2={Qn:.3f} nested Q2={Qg:.3f}')

RES=pd.DataFrame(rows); VIP=pd.DataFrame(vip_tab); TWO=pd.DataFrame(twostage)
fn=(_WORK + '/20260706_31MT_LitDeskriptoren_Analyse.xlsx')
with pd.ExcelWriter(fn,engine='openpyxl',mode='a',if_sheet_exists='replace') as w:
    RES.to_excel(w,sheet_name='PLS-DA_reduziert',index=False)
    VIP.to_excel(w,sheet_name='VIP_reduziert',index=False)
    TWO.to_excel(w,sheet_name='VIP_zweistufig',index=False)
RES.to_pickle((_WORK + '/lit_results.pkl')); TWO.to_pickle((_WORK + '/lit_twostage.pkl'))
print('\nGeschrieben.')
