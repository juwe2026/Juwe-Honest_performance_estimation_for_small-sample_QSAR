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
import sys, json; sys.path.insert(0,(_ROOT + '/code'))
import core
from sklearn.metrics import roc_auc_score
NLVCV = {'H.I.':4,'S.A.':4,'S.Pneu':3,'S.Pyo':1,'P.A.':4}
out={}
for tag,cols in (("D13",core.RED13),("D12",core.RED12)):
    out[tag]={}
    for s in core.STRAINS:
        n=NLVCV[s]; y=core.y_of(s); yh=core.loo(core.X_of(cols),y,n)
        c=core.confusion(cols,s,n)
        out[tag][s]=dict(nlv=n,q2=float(core.q2(y,yh)),auc=float(roc_auc_score(y,yh)),**c)
        print("%s %-7s nLV=%d Q2 %.3f AUC %.3f  %d/%d/%d/%d"%(tag,s,n,out[tag][s]["q2"],out[tag][s]["auc"],c["TP"],c["FN"],c["FP"],c["TN"]))
json.dump(out,open(_WORK + '/nlvcv.json','w'),indent=1)
