# Greift der Rueckfallzweig in den Gaucher-Armen (Tabelle S12, Additional file 13)?
# Gezaehlt mit instrumentierter Kopie von gaucher_pipe.py; die Konstante 0.5 steht
# bereits im korrigierten Code, der Zaehler aendert das Ergebnis nicht.
import sys, os, pickle, time
import numpy as np
sys.path.insert(0,'/home/claude/rev/auditcode_g')
import gaucher_pipe as gp
from sklearn.metrics import roc_auc_score
ROOT='/home/claude/rev/QSAR_selection_leakage'
d=pickle.load(open(ROOT+'/results/raw_permutations/gaucher_res.pkl','rb'))
groups=d['groups']; reps_u=d['reps_u']; y=d['y']
X,_,mz,_=gp.load(ROOT+'/data/Gaucherdata_paper39_median_normalised.csv')
Xc,_,_=gp.clean(X,mz)
print('Gaucher: %d Proben, %d Variablen nach Reinigung, %d aktiv'%(Xc.shape[0],Xc.shape[1],int(y.sum())),flush=True)
NP=int(os.environ.get('NP','100'))
for route in ('A1','A1u','A2'):
    gp.FB[0]=0; t=time.time()
    yh=gp.loo_nested(Xc,y,route,groups,reps_u); obs=gp.FB[0]
    oq=gp.q2(y,yh); oa=roc_auc_score(y,yh)
    pf=0; ph=0
    for k in range(NP):
        yp=np.random.default_rng([909,k]).permutation(y)
        gp.FB[0]=0
        gp.loo_nested(Xc,yp,route,groups,reps_u)
        pf+=gp.FB[0]; ph+= 1 if gp.FB[0] else 0
    print('%-4s Q2 %7.4f AUC %6.4f | echte Labels %d/%d Folds | permutiert %.3f %% der Folds, %.1f %% der Permutationen  (%.0fs)'
          %(route,oq,oa,obs,len(y),100.0*pf/(NP*len(y)),100.0*ph/NP,time.time()-t),flush=True)
