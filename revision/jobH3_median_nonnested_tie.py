# Median-Regel (Table S16): die vier Arme ohne genestete Auswahl mit Gleichstandsbehandlung neu
# (gleiche Seeds 9101-9104 wie jobH_modal.py); die genesteten Arme stehen schon korrigiert in H_median_fixed_*.json.
from common import *
import core, json
k=sys.argv[1]; NP=2000
MODAL={'HI':['Carvacrol','Citral','Citronellal','Menthol','Thymol','Thymoquinone','β-Thujaplicin'],
       'SPyo':['Carvacrol','Citral','Thymol','Thymoquinone','β-Thujaplicin']}
y=np.array([1.0 if n in MODAL[k] else 0.0 for n in core.COMPOUNDS]); assert y.sum()==len(MODAL[k])
X12=core.X_of(core.RED12); L25=sorted(set(core.LIT25)); X25=core.X_of(L25)
ARMS=[('RouteB_12_1LV',lambda yy: loo_pls(X12,yy,1),9101),('RouteB_12_3LV',lambda yy: loo_pls(X12,yy,3),9102),
      ('RouteB_25_1LV',lambda yy: loo_pls(X25,yy,1),9103),('noSelection_2208',lambda yy: loo_pls(X,yy,1),9104)]
p_out='/home/claude/rev/out/H_median_fixed_%s.json'%k
out=json.load(open(p_out))
for nm,f,seed in ARMS:
    old=out[nm]; t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    assert abs(oq-old['q2'])<1e-9 and abs(oa-old['auc'])<1e-9, (nm, oq, old['q2'])
    for j in range(NP):
        yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    pq,pa=pval(qs,oq),pval(as_,oa)
    out[nm]=dict(old,p_q2=pq,p_auc=pa,p_q2_strict=old['p_q2'],p_auc_strict=old['p_auc'],ties=True)
    print(k,nm,'Q2 p %.4f (strikt %.4f) | AUC p %.4f (strikt %.4f) %.0fs'%(pq,old['p_q2'],pa,old['p_auc'],time.time()-t),flush=True)
    json.dump(out,open(p_out,'w'),indent=1,default=float)
print('fertig',k,flush=True)
