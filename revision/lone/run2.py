# Zweiter MHK-Datensatz (Lone 2019). Wie run_lone.py, aber mit Cache pro Arm (neustartfest).
import sys, os, json, time
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import nested_generic as ng, altorder as ao, unsup_rep as ur
os.environ['QSAR_RECALC_DIR']='/home/claude/rev/lone/recalc'

RULE=sys.argv[1]            # 'A' oder 'C'
CACHE='/home/claude/rev/lone/cache'; os.makedirs(CACHE,exist_ok=True)
NPERM=2000; NPERM_RF=500; NPERM_HITS=1000
SEED={'A1':9201,'A1u':9202,'A2':9203,'nosel':9204,'nofilter':9205,'rf':9206,'hits':9207}
if RULE=='C': SEED={k:v+100 for k,v in SEED.items()}
THR={'A':{'EC':8,'SA':8,'PA':16,'AB':8},'C':{'EC':16,'SA':16,'PA':16,'AB':16}}[RULE]

rows=json.load(open('/home/claude/rev/lone/rows.json'))
X=pd.read_csv('/home/claude/rev/lone/desc_pool.csv'); assert list(X['No'].astype(str))==[r['No'] for r in rows]
DESC=[c for c in X.columns if c!='No']; X=X[DESC].values.astype(float); N=X.shape[0]
Y={}
for k in ['EC','SA','PA','AB']:
    Y[k]=np.array([1.0 if (float(r['MIC_'+k])<=THR[k] and r['censored_'+k]!='yes') else 0.0 for r in rows])
print('Regel',RULE,'Deskriptoren',X.shape,'Aktive',{k:int(v.sum()) for k,v in Y.items()},flush=True)

def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)
def pval(null,obs): null=np.asarray(null); return float((1+np.sum(null>=obs))/(1+len(null)))
def mcc(y,yh,thr=0.5):
    p=(yh>=thr).astype(float); tp=int(((p==1)&(y==1)).sum()); fn=int(((p==0)&(y==1)).sum())
    fp=int(((p==1)&(y==0)).sum()); tn=int(((p==0)&(y==0)).sum())
    den=np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return dict(TP=tp,FN=fn,FP=fp,TN=tn,mcc=float((tp*tn-fp*fn)/den) if den else 0.0)
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
t=time.time(); G,_=ao.build_groups(X); REPS=np.asarray(ur.centrality_reps(X,G))
print('Gruppen',len(G),'(%.0fs)'%(time.time()-t),flush=True)
R=rankdata(X,axis=0)

def cached(name,fn):
    p=os.path.join(CACHE,'%s_%s.json'%(RULE,name))
    if os.path.exists(p):
        print('cache',name,flush=True); return json.load(open(p))
    t=time.time(); v=fn(); v['seconds']=round(time.time()-t,1)
    json.dump(v,open(p,'w'),indent=1,default=float); return v

def perm(f,y,seed,npm=NPERM):
    yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    for j in range(npm):
        yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp)
        qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    return dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),nperm=npm,**mcc(y,yh))

out={'rule':RULE,'thresholds':THR,'n_desc':int(X.shape[1]),'n_groups':len(G),
     'n_active':{k:int(v.sum()) for k,v in Y.items()},'nperm':NPERM,'seeds':SEED}
for k,y in Y.items():
    D=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    out[k]={}
    arms=[('A1',lambda yy: ao.alt_nested_cv(D,G,yy,1)[2],SEED['A1']),
          ('A1u',lambda yy: ur.unsup_nested_cv(D,REPS,yy,1)[2],SEED['A1u']),
          ('A2',lambda yy: ng.nested_cv(D,yy,1)[2],SEED['A2']),
          ('noSelection_all',lambda yy: loo_pls(X,yy,1),SEED['nosel']),
          ('noFilter_reps',lambda yy: loo_pls(X[:,REPS],yy,1),SEED['nofilter'])]
    for nm,f,sd in arms:
        out[k][nm]=cached('%s_%s'%(k,nm),lambda f=f,sd=sd,y=y: perm(f,y,sd))
        r=out[k][nm]
        print('%s %-16s Q2 %6.3f p %.4f AUC %.3f p %.4f MCC %6.3f %ss'%(k,nm,r['q2'],r['p_q2'],r['auc'],r['p_auc'],r['mcc'],r.get('seconds')),flush=True)
        json.dump(out,open('/home/claude/rev/lone/res_%s.json'%RULE,'w'),indent=1,default=float)
    def rf_arm(y=y):
        def rf_oob(yy,seed=0):
            rf=RandomForestClassifier(300,min_samples_leaf=1,max_features='sqrt',oob_score=True,n_jobs=1,random_state=seed).fit(X,yy)
            return rf.oob_decision_function_[:,1]
        po=rf_oob(y); oq=q2(y,po); oa=roc_auc_score(y,po); qs=[];as_=[]
        for j in range(NPERM_RF):
            yp=np.random.default_rng([SEED['rf'],j]).permutation(y); p=rf_oob(yp); qs.append(q2(yp,p)); as_.append(roc_auc_score(yp,p))
        return dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),nperm=NPERM_RF,**mcc(y,po))
    out[k]['randomForest_oob']=cached('%s_rf'%k,rf_arm); r=out[k]['randomForest_oob']
    print('%s %-16s Q2 %6.3f p %.4f AUC %.3f p %.4f %ss'%(k,'randomForest',r['q2'],r['p_q2'],r['auc'],r['p_auc'],r.get('seconds')),flush=True)
    def hits_arm(y=y):
        def hits(yy):
            e,p=ng.eff_filter(R,yy,N); ps=(e>0.3)&(p<=0.05)
            a1=np.array([g[int(np.argmax(e[g]))] for g in G])
            return int(ps.sum()),int(ps[a1].sum()),int(ps[REPS].sum())
        ob=hits(y); nl=np.array([hits(np.random.default_rng([SEED['hits'],j]).permutation(y)) for j in range(NPERM_HITS)])
        d={}
        for i,nm in enumerate(['A2','A1','A1u']):
            ntests=X.shape[1] if nm=='A2' else len(G)
            m_=nl[:,i].mean(); s_=nl[:,i].std(ddof=1); q_=m_/ntests
            d[nm]=dict(n_tests=ntests,obs=ob[i],null_mean=float(m_),null_sd=float(s_),rate=float(q_),
                infl_nominal=float(s_/np.sqrt(ntests*.05*.95)),infl_own=float(s_/np.sqrt(ntests*q_*(1-q_))) if 0<q_<1 else None,
                p=pval(nl[:,i],ob[i]))
        return d
    out[k]['hits']=cached('%s_hits'%k,hits_arm)
    print(k,'hits',{n:(v['obs'],round(v['infl_own'],2) if v.get('infl_own') else None,round(v['p'],4)) for n,v in out[k]['hits'].items() if isinstance(v,dict)},flush=True)
    json.dump(out,open('/home/claude/rev/lone/res_%s.json'%RULE,'w'),indent=1,default=float)
print('fertig',RULE,flush=True)
