# Median rule (identical with the modal rule in all five strains, see
# the check in the revision log): the three nested Route A arms recomputed with the
# corrected fallback branch (0.5 instead of ytr.mean()), same seeds 9105-9107,
# 2,000 permutations. Route B and the arm without selection have no such branch
# and are taken from H_modal_<k>.json.
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
# These scripts imported core.py from the archive as submitted, whose default Route B file is the submitted one.
_osp.environ.setdefault('QSAR_ROUTEB_XLSX', _ROOT + '/data/RouteB_descriptor_data_submitted.xlsx')
import sys, os, json, time
os.environ.setdefault('QSAR_ROUTEA_XLSX',(_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
import core   # order of the compounds, from the original code (before the fixcode path)
assert _osp.path.realpath(_osp.path.dirname(core.__file__)) == _osp.path.realpath(_ROOT + '/code'), core.__file__
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,(_ROOT + '/code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur): assert (_ROOT + '/code') in m.__file__, m.__file__
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
k=sys.argv[1]; NP=2000
H=json.load(open((_OUT + '/H_modal_%s.json')%k))
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
p_out=(_OUT + '/H_median_fixed_%s.json')%k
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
    print(k,nm,'Q2 %.4f (old %.4f) p %.4f (old %.4f) | AUC p %.4f (old %.4f)  %.0fs'%(
        oq,old['q2'],out[nm]['p_q2'],old['p_q2'],out[nm]['p_auc'],old['p_auc'],time.time()-t),flush=True)
print('done',k,flush=True)
