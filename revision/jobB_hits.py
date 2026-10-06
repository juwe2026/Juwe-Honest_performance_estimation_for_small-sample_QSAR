# G2-6 / G2-9: tie-break sensitivity with own nulls; route-specific calibration of inflation factors
from common import *
import nested_generic as ng
H=np.load('/home/claude/rev/QSAR_selection_leakage/results/raw_permutations/hits_three_variants.npy',allow_pickle=True).item()
MAP={'HI':'H.I.','SA':'S.A.','SPneu':'S.Pneu','SPyo':'S.Pyo','PA':'P.A.'}
cal={}
for k in KEYS:
    h=H[MAP[k]]; cal[k]={}
    for key,n in [('A2',2208),('A1s',182),('A1u',182)]:
        nl=np.asarray(h['null_'+key],float); m=nl.mean(); sd=nl.std(ddof=1); q=m/n
        cal[k][key]=dict(obs=int(h['obs_'+key]),null_mean=m,null_sd=sd,rate=q,
            infl_nominal=sd/np.sqrt(n*.05*.95), infl_own=sd/np.sqrt(n*q*(1-q)), p=pval(nl,h['obs_'+key]))
for key in ['A2','A1s','A1u']:
    print(key,'rate %.3f-%.3f'%(min(cal[k][key]['rate'] for k in KEYS),max(cal[k][key]['rate'] for k in KEYS)),
      'nominal %.2f-%.2f'%(min(cal[k][key]['infl_nominal'] for k in KEYS),max(cal[k][key]['infl_nominal'] for k in KEYS)),
      'own %.2f-%.2f'%(min(cal[k][key]['infl_own'] for k in KEYS),max(cal[k][key]['infl_own'] for k in KEYS)))
    for k in KEYS: print('   ',k,'nom %.2f own %.2f'%(cal[k][key]['infl_nominal'],cal[k][key]['infl_own']))
# --- groupings: tie-break variants with their own null (1000 perms, seed [8202,j])
from scipy.stats import rankdata
R=rankdata(X,axis=0); p=X.shape[1]
Cp=np.abs(np.nan_to_num(np.corrcoef(X,rowvar=False))); Cs=np.abs(np.nan_to_num(np.corrcoef(R,rowvar=False)))
C=np.maximum(Cp,Cs); np.fill_diagonal(C,np.nan)
def groups(tb, rng=None):
    rem=np.ones(p,bool); deg=ADJ.sum(1).astype(int); Gs=[]; var=X.var(0)
    while rem.any():
        d=np.where(rem,deg,-1); mx=d.max()
        if mx<=0: Gs+=[[int(i)] for i in np.where(rem)[0]]; break
        cand=np.where(d==mx)[0]
        hub=cand[np.argmax(var[cand])] if tb=='var' else (cand[0] if tb=='first' else rng.choice(cand))
        mem=np.unique(np.concatenate([[hub],np.where(ADJ[hub]&rem)[0]])); Gs.append(mem.tolist())
        for m_ in mem: deg-=ADJ[m_].astype(int)
        rem[mem]=False
    return Gs
def creps(Gs):
    return np.array([g[0] if len(g)==1 else g[int(np.nanargmax(np.nanmean(C[np.ix_(g,g)],1)))] for g in Gs])
GV={'variance (as published)':groups('var'),'first index':groups('first')}
rng=np.random.default_rng(8203)
for r in range(100): GV['random_%03d'%r]=groups('rand',rng)
RP={g:creps(v) for g,v in GV.items()}
assert len(GV['variance (as published)'])==182
NP=1000
pass_obs={k:None for k in KEYS}; nulls={g:{k:[] for k in KEYS} for g in GV}
for k in KEYS:
    y=Y[k]; eff,pv=ng.eff_filter(R,y,N); pass_obs[k]=(eff>0.3)&(pv<=0.05)
    for j in range(NP):
        yp=np.random.default_rng([8202,j]).permutation(y); e,pv_=ng.eff_filter(R,yp,N); ps=(e>0.3)&(pv_<=0.05)
        for g in GV: nulls[g][k].append(int(ps[RP[g]].sum()))
res={}
for g in GV:
    res[g]={'n_groups':len(GV[g])}
    for k in KEYS:
        o=int(pass_obs[k][RP[g]].sum()); nl=np.array(nulls[g][k]); res[g][k]=dict(obs=o,null_mean=nl.mean(),null_sd=nl.std(ddof=1),p=pval(nl,o))
for g in ['variance (as published)','first index']:
    print(g,res[g]['n_groups'],[(k,res[g][k]['obs'],round(res[g][k]['p'],3)) for k in KEYS])
rand=[g for g in GV if g.startswith('random')]
summ={}
for k in KEYS:
    ps=np.array([res[g][k]['p'] for g in rand]); ob=np.array([res[g][k]['obs'] for g in rand])
    summ[k]=dict(obs_min=int(ob.min()),obs_max=int(ob.max()),obs_median=float(np.median(ob)),p_min=ps.min(),p_max=ps.max(),p_median=float(np.median(ps)),
                 frac_p_le_0125=float(np.mean(ps<=0.0125)),frac_p_le_05=float(np.mean(ps<=0.05)))
    print(k,summ[k])
ng_=[res[g]['n_groups'] for g in rand]; print('random n_groups',min(ng_),max(ng_),np.median(ng_))
dump('B_hits',dict(calibration=cal,tiebreak_fixed={g:res[g] for g in ['variance (as published)','first index']},random_summary=summ,random_ngroups=ng_,
     seeds=dict(null='default_rng([8202,j]), 1000 permutations',random_tiebreak='default_rng(8203), 100 groupings')))
# deposit groups
import csv
with open('/home/claude/rev/out/groups_182.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['group','size','hub_order','representative_A1u_centrality','member'])
    for gi,g in enumerate(GV['variance (as published)']):
        rep=RP['variance (as published)'][gi]
        for m_ in g: w.writerow([gi+1,len(g),gi+1,DESC[rep],DESC[m_]])
