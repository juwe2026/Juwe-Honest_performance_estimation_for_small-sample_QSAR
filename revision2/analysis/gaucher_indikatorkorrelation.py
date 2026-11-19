"""Post hoc, 05.11.2026: empirical correlation of the hit indicators under label permutation, Gaucher data set (code of Revision 1).
Filter as in Revision 1 (gaucher_pipe.eff_filter on ranks, |r_rb| > 0.3 and p <= 0.05), 3000 permutations, seed 1.
Output: null rate q, calibrated factor F, mean indicator correlation per |rho| class against the leading term k rho^2.
Usage: python -I gaucher_indikatorkorrelation.py <archive_root>   (archive_root: root folder of this archive)"""
import sys, os, numpy as np
from scipy.stats import rankdata, norm
sys.path.insert(0, os.path.join(sys.argv[1], 'code')); import gaucher_pipe as gp
X,y,mz,_=gp.load(os.path.join(sys.argv[1], 'data/Gaucherdata_paper39_median_normalised.csv')); X,_,_=gp.clean(X,mz)
n,p=X.shape; R=rankdata(X,axis=0)
C=np.corrcoef(R,rowvar=False); iu=np.triu_indices(p,1); rho=C[iu]
print('n',n,'p',p,'mean rho^2',np.mean(rho**2).round(4),'share of columns with ties',np.mean([len(np.unique(X[:,j]))<n for j in range(p)]).round(3))
print('distribution |rho|: quantiles', np.quantile(np.abs(rho),[.5,.9,.99]).round(3), 'share |rho|>0.7', np.mean(np.abs(rho)>0.7).round(4))
rng=np.random.default_rng(1); NP=3000; I=np.zeros((NP,p),bool)
for b in range(NP):
    yp=rng.permutation(y); eff,pv=gp.eff_filter(R,yp,n); I[b]=(eff>0.3)&(pv<=0.05)
q=I.mean(); cnt=I.sum(1); F=cnt.std(ddof=1)/np.sqrt(p*q*(1-q))
print('q',q.round(4),'F empirical',F.round(3))
Ic=I.astype(float)-I.mean(0); s=Ic.std(0); ok=s>0
Ci=np.corrcoef(I[:,ok].astype(float),rowvar=False); ius=np.triu_indices(ok.sum(),1)
rr=C[np.ix_(ok,ok)][ius]; ci=Ci[ius]
k=2*norm.ppf(0.975)**2*norm.pdf(norm.ppf(0.975))**2/(0.05*0.95)
for lo,hi in [(0,.1),(.1,.2),(.2,.3),(.3,.5),(.5,.7),(.7,.9),(.9,1.01)]:
    m=(np.abs(rr)>=lo)&(np.abs(rr)<hi)
    if m.sum(): print(f'|rho| {lo:.1f}-{hi:.1f}: pairs {m.sum():6d}  indicator correlation, empirical {ci[m].mean():.4f}  approximation k*rho^2 {np.mean(k*rr[m]**2):.4f}')
print('mean empirical',ci.mean().round(4),'approximation',np.mean(k*rr**2).round(4))
