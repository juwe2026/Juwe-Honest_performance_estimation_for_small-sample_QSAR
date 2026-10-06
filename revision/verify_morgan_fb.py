# Pruefung der Angaben "81 bis 88 % der Folds" und "mehr als neun von zehn Nullanpassungen" (Morgan-Arm, Seed 8602)
exec(open('/home/claude/rev/jobE2_morgan_fixed.py').read().split("E=json.load")[0])
import json
NPV=int(sys.argv[1]) if len(sys.argv)>1 else 300
res={}
def count_empty(Xm,y):
    c=0
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False; R=rankdata(Xm[tr],axis=0)
        eff,pv=ng.eff_filter(R,y[tr],30)
        if not np.any((eff>0.3)&(pv<=0.05)): c+=1
    return c
for k in KEYS:
    e=[count_empty(MB,np.random.default_rng([8602,j]).permutation(Y[k])) for j in range(NPV)]
    e=np.array(e); res[k]=dict(nperm=NPV,frac_folds=float(e.sum()/(31*NPV)),frac_perm_any=float((e>0).mean()),obs=count_empty(MB,Y[k]))
    print(k,res[k],flush=True)
json.dump(res,open('/home/claude/rev/out/morgan_fb_verify.json','w'),indent=1)
