# Zweiter MHK-Datensatz, Routen A1/A1u/A2 mit KORRIGIERTER Rueckfallvorhersage (Konstante 0.5
# statt Mittelwert der Trainingslabels). Die Arme ohne Filter und der Random Forest sind nicht
# betroffen und werden aus dem alten Cache uebernommen.
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/fixcode')
sys.path.insert(1,'/home/claude/rev/QSAR_selection_leakage/code')
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur):
    assert '/home/claude/rev/fixcode' in m.__file__, m.__file__

RULE=sys.argv[1]
CACHE='/home/claude/rev/lone/cache_fix'; os.makedirs(CACHE,exist_ok=True)
NPERM=2000
SEED={'A1':9201,'A1u':9202,'A2':9203}
if RULE=='C': SEED={k:v+100 for k,v in SEED.items()}
THR={'A':{'EC':8,'SA':8,'PA':16,'AB':8},'C':{'EC':16,'SA':16,'PA':16,'AB':16}}[RULE]

rows=json.load(open('/home/claude/rev/lone/rows.json'))
X=pd.read_csv('/home/claude/rev/lone/desc_pool.csv')
DESC=[c for c in X.columns if c!='No']; X=X[DESC].values.astype(float); N=X.shape[0]
Y={k:np.array([1.0 if (float(r['MIC_'+k])<=THR[k] and r['censored_'+k]!='yes') else 0.0 for r in rows])
   for k in ['EC','SA','PA','AB']}
print('Regel',RULE,'korrigierter Rueckfallwert','Aktive',{k:int(v.sum()) for k,v in Y.items()},flush=True)

def q2(y,yh): return 1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2)
exec(open('/home/claude/rev/pval_tie.py').read())   # Gleichstandsbehandlung, Revision 03.10.2026
def mcc(y,yh,thr=0.5):
    p=(yh>=thr).astype(float); tp=int(((p==1)&(y==1)).sum()); fn=int(((p==0)&(y==1)).sum())
    fp=int(((p==1)&(y==0)).sum()); tn=int(((p==0)&(y==0)).sum())
    den=np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return dict(TP=tp,FN=fn,FP=fp,TN=tn,mcc=float((tp*tn-fp*fn)/den) if den else 0.0)

FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
G,_=ao.build_groups(X); REPS=np.asarray(ur.centrality_reps(X,G))
print('Gruppen',len(G),flush=True)

out={'rule':RULE,'change':'empty-selection fallback 0.5 instead of mean of training labels',
     'n_groups':len(G),'n_desc':int(X.shape[1]),'seeds':SEED,'nperm':NPERM}
for k,y in Y.items():
    D=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    out[k]={}
    for nm,f in [('A1',lambda yy: ao.alt_nested_cv(D,G,yy,1)[2]),
                 ('A1u',lambda yy: ur.unsup_nested_cv(D,REPS,yy,1)[2]),
                 ('A2',lambda yy: ng.nested_cv(D,yy,1)[2])]:
        p=os.path.join(CACHE,'%s_%s_%s.json'%(RULE,k,nm))
        if os.path.exists(p):
            out[k][nm]=json.load(open(p)); print('cache',k,nm,flush=True)
        else:
            t=time.time(); sd=SEED[nm]
            yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
            for j in range(NPERM):
                yp=np.random.default_rng([sd,j]).permutation(y); h=f(yp)
                qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
            out[k][nm]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),nperm=NPERM,
                            seconds=round(time.time()-t,1),**mcc(y,yh))
            json.dump(out[k][nm],open(p,'w'),indent=1,default=float)
        r=out[k][nm]
        print('%s %-4s Q2 %6.3f p %.4f AUC %.3f p %.4f  %ss'%(k,nm,r['q2'],r['p_q2'],r['auc'],r['p_auc'],r.get('seconds')),flush=True)
        json.dump(out,open('/home/claude/rev/lone/res_%s_fix.json'%RULE,'w'),indent=1,default=float)
print('fertig',RULE,flush=True)
