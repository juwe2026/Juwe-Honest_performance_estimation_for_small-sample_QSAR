# Median rule (Table S16 of Revision 1): the four arms without nested selection recomputed with treatment of ties
# (same seeds 9101-9104 as jobH_modal.py), the nested arms are already corrected in H_median_fixed_*.json.
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
from common import *
import core, json
k=sys.argv[1]; NP=2000
MODAL={'HI':['Carvacrol','Citral','Citronellal','Menthol','Thymol','Thymoquinone','β-Thujaplicin'],
       'SPyo':['Carvacrol','Citral','Thymol','Thymoquinone','β-Thujaplicin']}
y=np.array([1.0 if n in MODAL[k] else 0.0 for n in core.COMPOUNDS]); assert y.sum()==len(MODAL[k])
X12=core.X_of(core.RED12); L25=sorted(set(core.LIT25)); X25=core.X_of(L25)
ARMS=[('RouteB_12_1LV',lambda yy: loo_pls(X12,yy,1),9101),('RouteB_12_3LV',lambda yy: loo_pls(X12,yy,3),9102),
      ('RouteB_25_1LV',lambda yy: loo_pls(X25,yy,1),9103),('noSelection_2208',lambda yy: loo_pls(X,yy,1),9104)]
p_out=(_OUT + '/H_median_fixed_%s.json')%k
out=json.load(open(p_out))
for nm,f,seed in ARMS:
    old=out[nm]; t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    assert abs(oq-old['q2'])<1e-9 and abs(oa-old['auc'])<1e-9, (nm, oq, old['q2'])
    for j in range(NP):
        yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    pq,pa=pval(qs,oq),pval(as_,oa)
    out[nm]=dict(old,p_q2=pq,p_auc=pa,p_q2_strict=old['p_q2'],p_auc_strict=old['p_auc'],ties=True)
    print(k,nm,'Q2 p %.4f (strict %.4f) | AUC p %.4f (strict %.4f) %.0fs'%(pq,old['p_q2'],pa,old['p_auc'],time.time()-t),flush=True)
    json.dump(out,open(p_out,'w'),indent=1,default=float)
print('done',k,flush=True)
