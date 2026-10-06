# KORRIGIERTE SEEDS: die erste Fassung hat fuer alle Arme 702 benutzt (das ist der Seed
# von Route B). Die veroeffentlichten Laeufe benutzen 59 (A1 top-15 nach Effektstaerke,
# prun_a1.py), 137 (A1 top-15 nach VIP mit einer Komponente, prun_vip1.py), 11 (A2
# top-15 nach Effektstaerke, prun15.py) und 23 (A2 top-15 nach VIP, prun_vip.py).
# Mit check_seeds.py gegen die ersten Nullwerte der abgelegten Dateien geprueft
# (Abweichung 0 fuer 31, 59, 137 und 11).
#
# Nachtrag zur Code-Korrektur: Die n:p-beschraenkten Auswahlen (Table 3) und die
# VIP-Auswahl (Table 4) laufen durch denselben fehlerhaften Rueckfallzweig wie die
# Routen A1/A1u/A2 und wurden bei der ersten Korrektur uebersehen.
# Hier: neu gerechnet mit der Konstanten 0.5 UND mit vollen 2.000 Permutationen
# (die abgelegten Dateien perm15_* und perm_vip_* haben nur 1.900 bis 1.982).
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,'/home/claude/rev/fixcode')
import nested_generic as ng, altorder as ao, unsup_rep as ur, a2vip
for m in (ng,ao,ur):
    assert '/home/claude/rev/fixcode' in m.__file__, m.__file__
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

NP=int(os.environ.get('NP','2000'))
ARMS=os.environ.get('ARMS','A1_top15,A1_viptop1,A1_top6,A2_top15,A2_viptop').split(',')
CACHE='/home/claude/rev/out/cache_topk_seeds'; os.makedirs(CACHE,exist_ok=True)
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
SEED={'A1_top15':59,'A1_viptop1':137,'A2_top15':11,'A2_viptop':23,'A1_viptop':73,'A1_top6':91}
print('Permutationen',NP,'| Arme',ARMS,flush=True)

def mcc(y,yh,thr=0.5):
    p=(yh>=thr).astype(float); tp=int(((p==1)&(y==1)).sum()); fn=int(((p==0)&(y==1)).sum())
    fp=int(((p==1)&(y==0)).sum()); tn=int(((p==0)&(y==0)).sum())
    den=np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return dict(TP=tp,FN=fn,FP=fp,TN=tn,mcc=float((tp*tn-fp*fn)/den) if den else 0.0)

for k in KEYS:
    y=Y[k]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    fns={'A1_top15' : lambda yy: ao.alt_nested_cv(D1,G,yy,1,15)[2],
         'A1_top6'  : lambda yy: ao.alt_nested_cv(D1,G,yy,1,6)[2],   # Revision: top-6 (prun_a1_k.py, Seed 91) mit korrigiertem Rueckfallzweig
         'A2_top15' : lambda yy: ng.nested_cv_a2_topk(D1,yy,1,15)[2],
         'A1_viptop': lambda yy: ao.alt_viptop(D1,G,yy,1,15,2)[2],
         'A1_viptop1': lambda yy: ao.alt_viptop(D1,G,yy,1,15,1)[2],
         'A2_viptop': lambda yy: a2vip.nested_cv_a2_viptop(D1,yy,1,15,2)[2]}
    for nm in ARMS:
        p=os.path.join(CACHE,'%s_%s.json'%(k,nm))
        if os.path.exists(p):
            print('cache',k,nm,flush=True); continue
        f=fns[nm]; t=time.time()
        yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
        for j in range(NP):
            yp=np.random.default_rng([SEED[nm],j]).permutation(y); h=f(yp)
            qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
        qs=np.array(qs); as_=np.array(as_)
        res=dict(q2=oq,auc=oa,
                 p_q2=pval(qs,oq),p_auc=pval(as_,oa),
                 nperm=NP,null_q2_median=float(np.median(qs)),null_auc_mean=float(as_.mean()),
                 blowups=int(np.sum(np.abs(qs)>10)),seconds=round(time.time()-t,1),**mcc(y,yh))
        json.dump(res,open(p,'w'),indent=1,default=float)
        print('%s %-10s Q2 %6.3f p %.4f AUC %.3f p %.4f  %.0fs'%(k,nm,oq,res['p_q2'],oa,res['p_auc'],time.time()-t),flush=True)
print('fertig',flush=True)
