# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 6 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import sys, os, time, numpy as np, pickle, gaucher_pipe as gp
from sklearn.metrics import roc_auc_score
route=sys.argv[1]; nlv=int(sys.argv[2]); target=int(sys.argv[3]); budget=float(sys.argv[4]) if len(sys.argv)>4 else 250.0
STORE=_os.path.join(_RES, f'gperm{nlv}_{route}.npz'); t0=time.time()
d=pickle.load(open(_os.path.join(_ROOT, 'results', 'raw_permutations', 'gaucher_res.pkl'),'rb')); groups=d['groups']; reps_u=d['reps_u']; y=d['y']
X,_,mz,_ = gp.load(_os.path.join(_DATA, 'Gaucherdata_paper39_median_normalised.csv')); Xc,_,_ = gp.clean(X,mz)
if os.path.exists(STORE):
    z=np.load(STORE); q2s=list(z['q2s']); aucs=list(z['aucs']); oq=float(z['obs_q2']); oa=float(z['obs_auc'])
else:
    yh=gp.loo_nested(Xc,y,route,groups,reps_u,nlv=nlv); oq=gp.q2(y,yh); oa=roc_auc_score(y,yh); q2s=[]; aucs=[]
def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
while len(q2s)<target and (time.time()-t0)<budget:
    k=len(q2s); yp=np.random.default_rng([424,k]).permutation(y)
    yh=gp.loo_nested(Xc,yp,route,groups,reps_u,nlv=nlv)
    q2s.append(gp.q2(yp,yh)); aucs.append(roc_auc_score(yp,yh))
    if len(q2s)%100==0: save()
save(); n=len(q2s)
pq=(1+np.sum(np.array(q2s)>=oq - 1e-9))/(1+n); pa=(1+np.sum(np.array(aucs)>=oa - 1e-9))/(1+n)
print(f'{route} nLV={nlv}: n={n} | Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pq:.4f} p(AUC)={pa:.4f} | {time.time()-t0:.0f}s')
