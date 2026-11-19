# Revision 03.10.2026: permutation p-values count ties (tolerance 1e-9), see CHANGELOG_revision.md
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
import sys, os, time, numpy as np, nested_generic as ng
short=sys.argv[1]; sheet=sys.argv[2]; target=int(sys.argv[3]); code=sys.argv[4]
budget=float(sys.argv[5]) if len(sys.argv)>5 else 255.0
STORE=f'{_WORK}/perm15_{code}.npz'; t0=time.time()
D=ng.build(sheet)
if os.path.exists(STORE):
    d=np.load(STORE); q2s=list(d['q2s']); aucs=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
else:
    oq,oa,_=ng.nested_cv_topk(D,D['Y'],1,top_k=15); q2s=[]; aucs=[]
def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
while len(q2s)<target and (time.time()-t0)<budget:
    k=len(q2s); yv=np.random.default_rng([11,k]).permutation(D['Y'])
    q,a,_=ng.nested_cv_topk(D,yv,1,top_k=15); q2s.append(q); aucs.append(a)
    if len(q2s)%50==0: save()
save(); n=len(q2s)
pQ=(1+np.sum(np.array(q2s)>=oq - 1e-9))/(1+n); pA=(1+np.sum(np.array(aucs)>=oa - 1e-9))/(1+n)
print(f'{short} top15: {n}/{target} | obs Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | {time.time()-t0:.0f}s')
