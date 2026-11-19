# Does the fallback branch apply in the Gaucher arms (Table S12 and Additional file 13 of Revision 1)?
# Counted with an instrumented copy of gaucher_pipe.py. The constant 0.5 is
# already in the corrected code, the counter does not change the result.
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)
import sys, os, pickle, time
import numpy as np
sys.path.insert(0,(_REV + '/auditcode_g'))
import gaucher_pipe as gp
from sklearn.metrics import roc_auc_score
ROOT=_ROOT
d=pickle.load(open(ROOT+'/results/raw_permutations/gaucher_res.pkl','rb'))
groups=d['groups']; reps_u=d['reps_u']; y=d['y']
X,_,mz,_=gp.load(ROOT+'/data/Gaucherdata_paper39_median_normalised.csv')
Xc,_,_=gp.clean(X,mz)
print('Gaucher: %d samples, %d variables after cleaning, %d active'%(Xc.shape[0],Xc.shape[1],int(y.sum())),flush=True)
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
    print('%-4s Q2 %7.4f AUC %6.4f | true labels %d/%d folds | permuted %.3f %% of the folds, %.1f %% of the permutations  (%.0fs)'
          %(route,oq,oa,obs,len(y),100.0*pf/(NP*len(y)),100.0*ph/NP,time.time()-t),flush=True)
