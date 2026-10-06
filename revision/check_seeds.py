# Welcher Seed steckt hinter welcher abgelegten Permutationsdatei?
# Pruefung: die ersten Permutationen mit dem ORIGINALCODE (fehlerhafter Rueckfallzweig)
# und dem vermuteten Seed nachrechnen und gegen die abgelegten Nullwerte stellen.
# Trifft es, stimmt die Seed-Zuordnung; dann muessen die korrigierten Laeufe denselben
# Seed benutzen, sonst mischt sich Monte-Carlo-Rauschen in die Differenz.
import sys, os, numpy as np
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
import nested_generic as ng, altorder as ao, unsup_rep as ur
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
assert 'QSAR_selection_leakage/code' in ng.__file__, ng.__file__

FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))

RP='/home/claude/rev/QSAR_selection_leakage/results/raw_permutations/'
ARMS={
 'perm_alt_%s.npz'  : (31, lambda D1,G,yy: ao.alt_nested_cv(D1,G,yy,1)[2]),
 'perm_uns_%s.npz'  : (47, None),
 'perm_%s.npz'      : (7,  lambda D1,G,yy: ng.nested_cv_a2(D1,yy,1)[2] if hasattr(ng,'nested_cv_a2') else None),
 'perm_a1e_%s.npz'  : (59, lambda D1,G,yy: ao.alt_nested_cv(D1,G,yy,1,15)[2]),
 'perm_a1vip1_%s.npz':(137,lambda D1,G,yy: ao.alt_viptop(D1,G,yy,1,15,1)[2]),
 'perm15_%s.npz'    : (11, lambda D1,G,yy: ng.nested_cv_a2_topk(D1,yy,1,15)[2]),
}
K='HI'; y=Y[K]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
for pat,(seed,fn) in ARMS.items():
    if fn is None: continue
    f=RP+pat%K
    if not os.path.exists(f):
        print('%-22s Datei fehlt'%pat); continue
    z=np.load(f,allow_pickle=True); q2s=np.asarray(z['q2s'],float)
    ok=[]
    for j in range(3):
        yp=np.random.default_rng([seed,j]).permutation(y)
        h=fn(D1,G,yp)
        if h is None: ok=None; break
        ok.append(abs(q2(yp,h)-q2s[j]))
    if ok is None:
        print('%-22s Funktion nicht verfuegbar'%pat); continue
    print('%-22s Seed %4d | Abweichung der ersten drei Nullwerte: %s'
          %(pat,seed,', '.join('%.2e'%d for d in ok)))
