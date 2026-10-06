import sys,json; sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
import numpy as np, gaucher_pipe as gp
from scipy.stats import rankdata
X,y,mz,_=gp.load('/home/claude/rev/QSAR_selection_leakage/data/Gaucherdata_paper39_median_normalised.csv'); X,mz,_=gp.clean(X,mz)
n,p=X.shape; G,_=gp.build_groups(X); reps=np.array(gp.centrality_reps(X,G))
R=rankdata(X,axis=0)
def counts(yy):
    e,pv=gp.eff_filter(R,yy,n); ps=(e>0.3)&(pv<=0.05)
    a1=np.array([g[int(np.argmax(e[g]))] for g in G])
    return int(ps.sum()), int(ps[a1].sum()), int(ps[reps].sum())
obs=counts(y); nul=np.array([counts(np.random.default_rng([8801,j]).permutation(y)) for j in range(1000)])
out={'n_var':p,'n_groups':len(G)}
for i,(nm,nt) in enumerate([('A2',p),('A1',len(G)),('A1u',len(G))]):
    nl=nul[:,i]; m=nl.mean(); sd=nl.std(ddof=1); q=m/nt
    out[nm]=dict(obs=obs[i],null_mean=m,null_sd=sd,rate=q,infl_nominal=sd/np.sqrt(nt*.05*.95),infl_own=sd/np.sqrt(nt*q*(1-q)),p=(1+np.sum(nl>=obs[i]))/1001)
    print(nm,{k:round(v,3) if isinstance(v,float) else v for k,v in out[nm].items()})
print(out['n_var'],out['n_groups'])
json.dump(out,open('/home/claude/rev/out/G_gaucher_hits.json','w'),default=float,indent=1)
