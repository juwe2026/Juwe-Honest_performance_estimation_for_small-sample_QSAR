# Check of the statements "81 to 88 % of the folds" and "more than nine of ten null fits" (Morgan arm, seed 8602)
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
exec(open((_REV + '/jobE2_morgan_fixed.py')).read().split("E=json.load")[0])
import json
NPV=int(sys.argv[1]) if len(sys.argv)>1 else 300
res={}
def count_empty(Xm,y):
    c=0
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False; R=rankdata(Xm[tr],axis=0)
        eff,pv=ng.eff_filter(R,y[tr],30)
        if not np.any((eff>0.3)&(pv<=0.05)): c+=1
    return c
for k in KEYS:
    e=[count_empty(MB,np.random.default_rng([8602,j]).permutation(Y[k])) for j in range(NPV)]
    e=np.array(e); res[k]=dict(nperm=NPV,frac_folds=float(e.sum()/(31*NPV)),frac_perm_any=float((e>0).mean()),obs=count_empty(MB,Y[k]))
    print(k,res[k],flush=True)
json.dump(res,open((_OUT + '/morgan_fb_verify.json'),'w'),indent=1)
