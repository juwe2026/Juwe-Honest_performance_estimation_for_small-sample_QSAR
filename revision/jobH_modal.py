# Sensitivity to the label definition: modal liquid-phase MIC (Werle 2026 BMCL Table 1, broth MIC <= 1024 = active)
from common import *
import core, altorder as ao, unsup_rep as ur, nested_generic as ng
from scipy.stats import rankdata
k=sys.argv[1]; NP=2000
MODAL={'HI':['Carvacrol','Citral','Citronellal','Menthol','Thymol','Thymoquinone','β-Thujaplicin'],
       'SPyo':['Carvacrol','Citral','Thymol','Thymoquinone','β-Thujaplicin']}
names=core.COMPOUNDS
y=np.array([1.0 if n in MODAL[k] else 0.0 for n in names]); assert y.sum()==len(MODAL[k])
d=D[k]
X12=core.X_of(core.RED12); L25=sorted(set(core.LIT25)); X25=core.X_of(L25)
ARMS=[('RouteB_12_1LV',lambda yy: loo_pls(X12,yy,1),9101),('RouteB_12_3LV',lambda yy: loo_pls(X12,yy,3),9102),
      ('RouteB_25_1LV',lambda yy: loo_pls(X25,yy,1),9103),('noSelection_2208',lambda yy: loo_pls(X,yy,1),9104),
      ('A1_nested',lambda yy: ao.alt_nested_cv(d,G,yy,1)[2],9105),('A1u_nested',lambda yy: ur.unsup_nested_cv(d,REPS_U,yy,1)[2],9106),
      ('A2_nested',lambda yy: ng.nested_cv(d,yy,1)[2],9107)]
out={'labels':'modal broth MIC <= 1024 (BMCL Table 1)','n_active':int(y.sum()),'actives':MODAL[k]}
R=rankdata(X,axis=0)
def hits(yy):
    e,p=ng.eff_filter(R,yy,N); ps=(e>0.3)&(p<=0.05)
    a1=np.array([g[int(np.argmax(e[g]))] for g in G]); return int(ps.sum()),int(ps[a1].sum()),int(ps[REPS_U].sum())
ob=hits(y); nl=np.array([hits(np.random.default_rng([9108,j]).permutation(y)) for j in range(1000)])
out['hits']={nm:dict(obs=ob[i],null_mean=float(nl[:,i].mean()),null_sd=float(nl[:,i].std(ddof=1)),p=pval(nl[:,i],ob[i])) for i,nm in enumerate(['A2','A1','A1u'])}
print(k,'hits',out['hits'],flush=True)
for nm,f,seed in ARMS:
    t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    for j in range(NP):
        yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    out[nm]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),**mcc_conf(y,yh))
    print(k,nm,'Q2 %.3f p %.4f AUC %.3f p %.4f MCC %.3f %.0fs'%(oq,out[nm]['p_q2'],oa,out[nm]['p_auc'],out[nm]['mcc'],time.time()-t),flush=True)
    dump('H_modal_'+k,out)
