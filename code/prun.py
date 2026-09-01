import sys, os, time, numpy as np, nested_generic as ng
short=sys.argv[1]; sheet=sys.argv[2]; target=int(sys.argv[3]); code=sys.argv[4]
budget=float(sys.argv[5]) if len(sys.argv)>5 else 240.0
STORE=f'/home/claude/perm_{code}.npz'; t0=time.time()
D=ng.build(sheet)
if os.path.exists(STORE):
    d=np.load(STORE); q2s=list(d['q2s']); aucs=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
else:
    oq,oa,_=ng.nested_cv(D,D['Y'],1); q2s=[]; aucs=[]
def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
while len(q2s)<target and (time.time()-t0)<budget:
    k=len(q2s)
    yv=np.random.default_rng([7,k]).permutation(D['Y'])
    q,a,_=ng.nested_cv(D,yv,1); q2s.append(q); aucs.append(a)
    if len(q2s)%50==0: save()
save(); n=len(q2s)
pQ=(1+np.sum(np.array(q2s)>=oq))/(1+n); pA=(1+np.sum(np.array(aucs)>=oa))/(1+n)
print(f'{short}: perms={n}/{target} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | obs Q2={oq:.3f} AUC={oa:.3f} | null-med Q2={np.median(q2s):.3f} AUC={np.median(aucs):.3f} | {time.time()-t0:.0f}s')
