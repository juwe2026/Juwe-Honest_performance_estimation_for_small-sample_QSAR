# CORRECTED VERSION of the Morgan arm with nested filter (jobE_misc.py, part 3).
# The own filter function set ytr.mean() when the selection was empty. With the Morgan bits
# no bit passes the filter for H. influenzae and S. pyogenes in most folds,
# so the branch already applies with the true labels and distorts the point estimate
# (AUC 0.000 = perfectly reversed ranking). Now the neutral constant 0.5,
# same seed 8602, 2,000 permutations. The arm without selection has no such
# branch and stays unchanged.
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
import sys, os
os.environ.setdefault('QSAR_ROUTEA_XLSX',(_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
import core
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,(_ROOT + '/code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur): assert (_ROOT + '/code') in m.__file__, m.__file__
from scipy.stats import rankdata
import json
NP=2000
MB=np.load((_OUT + '/morgan2_bits.npy'))
FB=[0]
def nested_filter(Xm,y):
    yh=np.zeros(N)
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False; ytr=y[tr]; R=rankdata(Xm[tr],axis=0)
        eff,pv=ng.eff_filter(R,ytr,30); keep=np.where((eff>0.3)&(pv<=0.05))[0]
        if len(keep)==0: yh[i]=0.5; FB[0]+=1; continue
        Xt=Xm[tr][:,keep]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1); ok=sd>1e-6; Xt=Xt[:,ok]; mu=mu[ok]; sd=sd[ok]
        if Xt.shape[1]==0: yh[i]=0.5; FB[0]+=1; continue
        m=PLSRegression(1,scale=False).fit((Xt-mu)/sd,ytr); yh[i]=m.predict((Xm[i:i+1][:,keep][:,ok]-mu)/sd).ravel()[0]
    return yh
E=json.load(open((_OUT + '/E_misc.json')))
out={'n_bits':int(MB.shape[1])}
for k in KEYS:
    y=Y[k]; FB[0]=0
    yh=nested_filter(MB,y); obs_fb=FB[0]; oq=q2(y,yh)
    try: oa=roc_auc_score(y,yh)
    except Exception: oa=float('nan')
    qs=[];as_=[]
    for j in range(NP):
        yp=np.random.default_rng([8602,j]).permutation(y); h=nested_filter(MB,yp)
        qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    old=E['morgan2'][k]['nested_filter']
    out[k]=dict(no_selection=E['morgan2'][k]['no_selection'],
                nested_filter=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),**mcc_conf(y,yh),
                                   empty_folds_observed=obs_fb,q2_old=old['q2'],auc_old=old['auc'],
                                   p_q2_old=old['p_q2'],p_auc_old=old['p_auc']))
    r=out[k]['nested_filter']
    print('%-6s empty folds %2d/31 | Q2 %.3f (old %.3f) p %.4f | AUC %.3f (old %.3f) p %.4f'%(
        k,obs_fb,oq,old['q2'],r['p_q2'],oa,old['auc'],r['p_auc']),flush=True)
    json.dump(out,open((_OUT + '/morgan2_fixed.json'),'w'),indent=1,default=float)
print('done',flush=True)
