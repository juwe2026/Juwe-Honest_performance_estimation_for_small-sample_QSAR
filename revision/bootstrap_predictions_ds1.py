# Observed predictions of every monoterpenoid model whose honest Q2 is compared in the manuscript (input of paired_bootstrap.py).
# No permutations, deterministic (random forest: random_state 0, out-of-bag). Same computations as jobE_misc.py, jobJ2_topk.py
# and jobD_rf.py; the Q2 of every model reproduces the published value. Paths as in the other revision scripts (README_revision.md).
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
import sys, os, json
os.environ.setdefault('QSAR_ROUTEA_XLSX',(_ROOT + '/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx'))
sys.path.insert(0,(_ROOT + '/code'))
exec(open((_REV + '/common.py')).read())
import core
for _m in ['nested_generic','altorder','unsup_rep']: sys.modules.pop(_m,None)
sys.path.insert(0,(_ROOT + '/code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur, a2vip
print('core', core.__file__, 'RED12', core.RED12)
from scipy.stats import rankdata
from sklearn.ensemble import RandomForestClassifier
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
S={'HI':'H.I.','SA':'S.A.','SPneu':'S.Pneu','SPyo':'S.Pyo','PA':'P.A.'}
def rf_oob(Xm,y):
    rf=RandomForestClassifier(300,min_samples_leaf=1,max_features='sqrt',oob_score=True,n_jobs=1,random_state=0).fit(Xm,y)
    return rf.oob_decision_function_[:,1]
XB=core.X_of(core.RED12); RED11=[c for c in core.RED12 if c!='Phenol_flag']; XB11=core.X_of(RED11)
L25=sorted(set(core.LIT25)); X25=core.X_of(L25)
out={}
for k in KEYS:
    y=Y[k]; assert np.array_equal(y, core.y_of(S[k])), k
    D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    P={}
    P['A1']=ao.alt_nested_cv(D1,G,y,1)[2]
    P['A1u']=ur.unsup_nested_cv(D1,REPS_U,y,1)[2]
    P['A2']=ng.nested_cv(D1,y,1)[2]
    P['NoSel2208']=loo_pls(X,y,1)
    P['Reps182']=loo_pls(X[:,REPS_U],y,1)
    P['RF2208']=rf_oob(X,y)
    P['B12']=core.loo(XB,y,1)
    P['B12_3LV']=core.loo(XB,y,3)
    P['B11_noPhenol']=core.loo(XB11,y,1)
    P['B25']=core.loo(X25,y,1)
    P['A1_top15']=ao.alt_nested_cv(D1,G,y,1,15)[2]
    P['A1_vip15']=ao.alt_viptop(D1,G,y,1,15,1)[2]
    P['A1_top6']=ao.alt_nested_cv(D1,G,y,1,6)[2]
    P['A2_top15']=ng.nested_cv_a2_topk(D1,y,1,15)[2]
    P['A2_vip15']=a2vip.nested_cv_a2_viptop(D1,y,1,15,2)[2]
    out[k]=dict(y=y.tolist(),pred={n:np.asarray(v).tolist() for n,v in P.items()})
    print(k,{n:round(q2(y,np.asarray(v)),3) for n,v in P.items()},flush=True)
json.dump(out,open((_OUT + '/bootstrap_predictions_ds1.json'),'w'))
