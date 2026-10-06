# Median-Regel (identisch mit der Modalwert-Regel in allen fuenf Staemmen, siehe
# Pruefung im Revisionsprotokoll): die drei genesteten Route-A-Arme neu mit dem
# korrigierten Rueckfallzweig (0.5 statt ytr.mean()), gleiche Seeds 9105-9107,
# 2.000 Permutationen. Route B und der Arm ohne Auswahl haben keinen solchen Zweig
# und werden aus H_modal_<k>.json uebernommen.
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
import core   # Reihenfolge der Verbindungen, aus dem Originalcode (vor dem fixcode-Pfad)
assert 'QSAR_selection_leakage/code' in core.__file__, core.__file__
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,'/home/claude/rev/fixcode')
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur): assert '/home/claude/rev/fixcode' in m.__file__, m.__file__
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
k=sys.argv[1]; NP=2000
H=json.load(open('/home/claude/rev/out/H_modal_%s.json'%k))
names=core.COMPOUNDS
y=np.array([1.0 if n in H['actives'] else 0.0 for n in names]); assert y.sum()==H['n_active']
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
d=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
ARMS=[('A1_nested',lambda yy: ao.alt_nested_cv(d,G,yy,1)[2],9105),
      ('A1u_nested',lambda yy: ur.unsup_nested_cv(d,REPS_U,yy,1)[2],9106),
      ('A2_nested',lambda yy: ng.nested_cv(d,yy,1)[2],9107)]
out=dict(H); out['labels']='median liquid-phase MIC <= 1024 (identical to the modal value for all five strains)'
p_out='/home/claude/rev/out/H_median_fixed_%s.json'%k
if os.path.exists(p_out): out=json.load(open(p_out))
for nm,f,seed in ARMS:
    if out.get(nm,{}).get('fixed'): print('cache',k,nm,flush=True); continue
    t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    for j in range(NP):
        yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    old=H[nm]
    out[nm]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),**mcc_conf(y,yh),fixed=True,
                 p_q2_old=old['p_q2'],p_auc_old=old['p_auc'],q2_old=old['q2'])
    json.dump(out,open(p_out,'w'),indent=1,default=float)
    print(k,nm,'Q2 %.4f (alt %.4f) p %.4f (alt %.4f) | AUC p %.4f (alt %.4f)  %.0fs'%(
        oq,old['q2'],out[nm]['p_q2'],old['p_q2'],out[nm]['p_auc'],old['p_auc'],time.time()-t),flush=True)
print('fertig',k,flush=True)
