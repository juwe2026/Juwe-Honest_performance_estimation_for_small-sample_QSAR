# Pruefung der Rekonstruktion: reproduziert nested_cv_a2_viptop die abgelegten
# Beobachtungswerte von perm_vip_*.npz? Der Rueckfallzweig ist in fixcode korrigiert,
# er greift bei den echten Labels aber praktisch nie, daher muessen die
# Beobachtungswerte uebereinstimmen, wenn die Rekonstruktion richtig ist.
import sys, os, glob
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,'/home/claude/rev/fixcode')
import nested_generic as ng, a2vip
from scipy.stats import rankdata
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
RP='/home/claude/rev/QSAR_selection_leakage/results/raw_permutations/'
for nlv_vip in (2,1):
    print('--- nlv_vip =',nlv_vip)
    for k in KEYS:
        y=Y[k]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
        q,a,_=a2vip.nested_cv_a2_viptop(D1,y,1,15,nlv_vip)
        z=np.load(RP+'perm_vip_%s.npz'%k,allow_pickle=True)
        oq=float(z['obs_q2']); oa=float(z['obs_auc'])
        print('%-6s Q2 %7.4f (Archiv %7.4f, d=%7.4f) | AUC %6.4f (Archiv %6.4f, d=%6.4f)'
              %(k,q,oq,q-oq,a,oa,a-oa))
