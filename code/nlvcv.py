import sys, json; sys.path.insert(0,'/home/claude/work/recalc')
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
json.dump(out,open("recalc/nlvcv.json","w"),indent=1)
