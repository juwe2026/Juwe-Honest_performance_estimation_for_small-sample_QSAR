import os, json
# G1-2 / G2-12: no-filter baselines, 1 LV PLS-DA, LOO, 2000 permutations, seed [8101,k]
from common import *
import nested_generic as ng
NPERM=int(sys.argv[1]) if len(sys.argv)>1 else 2000
R_all=None
def a1reps_nofilter(y):
    # Route A1 representatives (largest |effect| per group) recomputed in every fold, NO filter
    yh=np.zeros(N)
    for F in D['HI']['FOLDS']:
        tr=F['tr']; ytr=y[tr]; eff,pv=ng.eff_filter(F['R'],ytr,int(tr.sum()))
        reps=np.array([g[int(np.argmax(eff[g]))] for g in G])
        Xt=F['Xtr'][:,reps]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1); ok=sd>1e-6
        Xt=Xt[:,ok]; mu=mu[ok]; sd=sd[ok]
        m=PLSRegression(1,scale=False).fit((Xt-mu)/sd,ytr)
        yh[F['i']]=m.predict((F['Xte'][:,reps][:,ok]-mu)/sd).ravel()[0]
    return yh
ARMS={'all2208':lambda y: loo_pls(X,y), 'reps182_centrality':lambda y: loo_pls(X[:,REPS_U],y),
      'reps182_effectsize_nested':a1reps_nofilter}
out={}
_pp='/home/claude/rev/out/A_nofilter.json'
if os.environ.get('RESUME')=='1' and os.path.exists(_pp): out=json.load(open(_pp))   # Neustart nach Container-Wechsel
for arm,f in ARMS.items():
    out.setdefault(arm,{})
    for k in KEYS:
        if k in out[arm]: print('cache',arm,k,flush=True); continue
        t=time.time(); y=Y[k]; yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
        for j in range(NPERM):
            yp=np.random.default_rng([8101,j]).permutation(y); yhp=f(yp); qs.append(q2(yp,yhp)); as_.append(roc_auc_score(yp,yhp))
        c=mcc_conf(y,yh)
        out[arm][k]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),nperm=NPERM,seed='default_rng([8101,k])',**c)
        print(arm,k,'Q2 %.3f p %.4f AUC %.3f p %.4f MCC %.3f  %.0fs'%(oq,out[arm][k]['p_q2'],oa,out[arm][k]['p_auc'],c['mcc'],time.time()-t),flush=True)
        dump('A_nofilter',out)
