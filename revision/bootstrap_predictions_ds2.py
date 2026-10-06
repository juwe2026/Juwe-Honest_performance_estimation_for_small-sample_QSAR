# Observed predictions of the octapeptide models (main activity definition) compared in the manuscript (input of
# paired_bootstrap.py). Same computations as lone/run2_tie.py and lone/run3.py, no permutations, deterministic.
import sys, os, json
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
os.environ['QSAR_RECALC_DIR']='/home/claude/rev/lone/recalc'
sys.path.insert(0,'/home/claude/rev/fixcode'); sys.path.insert(1,'/home/claude/rev/QSAR_selection_leakage/code')
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestClassifier
import nested_generic as ng, altorder as ao, unsup_rep as ur
THR={'EC':8,'SA':8,'PA':16,'AB':8}
rows=json.load(open('/home/claude/rev/lone/rows.json'))
X=pd.read_csv('/home/claude/rev/lone/desc_pool.csv'); DESC=[c for c in X.columns if c!='No']; X=X[DESC].values.astype(float); N=X.shape[0]
Y={k:np.array([1.0 if (float(r['MIC_'+k])<=THR[k] and r['censored_'+k]!='yes') else 0.0 for r in rows]) for k in THR}
def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)
def loo_pls(Xm,y,nlv=1):
    yh=np.zeros(N)
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False
        Xt=Xm[tr]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1); ok=sd>1e-6
        Xt=Xt[:,ok]; mu=mu[ok]; sd=sd[ok]
        m=PLSRegression(n_components=min(nlv,Xt.shape[1],N-2),scale=False).fit((Xt-mu)/sd,y[tr])
        yh[i]=m.predict((Xm[i:i+1][:,ok]-mu)/sd).ravel()[0]
    return yh
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
G,_=ao.build_groups(X); REPS=np.asarray(ur.centrality_reps(X,G)); print('Gruppen',len(G))
out={}
for k,y in Y.items():
    D=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    rf=RandomForestClassifier(300,min_samples_leaf=1,max_features='sqrt',oob_score=True,n_jobs=1,random_state=0).fit(X,y)
    P={'A1':ao.alt_nested_cv(D,G,y,1)[2],'A1u':ur.unsup_nested_cv(D,REPS,y,1)[2],'A2':ng.nested_cv(D,y,1)[2],
       'NoSel1275':loo_pls(X,y,1),'Reps50':loo_pls(X[:,REPS],y,1),'RF1275':rf.oob_decision_function_[:,1]}
    out[k]=dict(y=y.tolist(),pred={n:np.asarray(v).tolist() for n,v in P.items()})
    print(k,{n:round(q2(y,np.asarray(v)),3) for n,v in P.items()},flush=True)
json.dump(out,open('/home/claude/rev/out/bootstrap_predictions_ds2.json','w'))
