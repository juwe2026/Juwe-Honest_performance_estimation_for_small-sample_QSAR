# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
import sys, os, time, numpy as np, nested_generic as ng, altorder as ao
short=sys.argv[1]; sheet=sys.argv[2]; target=int(sys.argv[3]); code=sys.argv[4]; topk=int(sys.argv[5])
budget=float(sys.argv[6]) if len(sys.argv)>6 else 250.0
STORE=f'/home/claude/perm_a1k{topk}_{code}.npz'; t0=time.time()
D=ng.build(sheet); groups=ao.get_groups(D['X'])
if os.path.exists(STORE):
    d=np.load(STORE); q2s=list(d['q2s']); aucs=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
else:
    oq,oa,_=ao.alt_nested_cv(D,groups,D['Y'],1,top_k=topk); q2s=[]; aucs=[]
def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
while len(q2s)<target and (time.time()-t0)<budget:
    k=len(q2s); yv=np.random.default_rng([91,k]).permutation(D['Y'])
    q,a,_=ao.alt_nested_cv(D,groups,yv,1,top_k=topk); q2s.append(q); aucs.append(a)
    if len(q2s)%200==0: save()
save(); n=len(q2s)
pQ=(1+np.sum(np.array(q2s)>=oq - 1e-9))/(1+n); pA=(1+np.sum(np.array(aucs)>=oa - 1e-9))/(1+n)
print(f'{short} A1-top{topk}: n={n} | Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | {time.time()-t0:.0f}s')
