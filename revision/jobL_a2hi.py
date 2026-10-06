# Route A2, H. influenzae: eigener Permutationslauf im abgelegten Stand.
# Die vier anderen Staemme stammen aus prun.py (Seed 7), H. influenzae aus
# perm_runner.py (Seed 42, ueber nested_h.py). Mit check_seed_a2c.py gegen den
# ersten Nullwert von perm_HI.npz geprueft: Seed 42 trifft exakt, 7 nicht.
# Hier mit korrigiertem Rueckfallzweig und 2.000 Permutationen nachgerechnet.
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,'/home/claude/rev/fixcode')
import nested_generic as ng
assert '/home/claude/rev/fixcode' in ng.__file__
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
NP=2000; SEED=42
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
y=Y['HI']; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
f=lambda yy: ng.nested_cv(D1,yy,1)[2]
t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
for j in range(NP):
    yp=np.random.default_rng([SEED,j]).permutation(y); h=f(yp)
    qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
qs=np.array(qs); as_=np.array(as_)
res=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),
         nperm=NP,seed=SEED,seconds=round(time.time()-t,1))
json.dump(res,open('/home/claude/rev/out/A2_HI_seed42.json','w'),indent=1,default=float)
print('A2 HI Seed 42: Q2 %.4f p %.4f | AUC %.4f p %.4f  (%.0fs)'%(oq,res['p_q2'],oa,res['p_auc'],res['seconds']),flush=True)
