from common import *
import nested_generic as ng, altorder as ao, unsup_rep as ur, core
from scipy.stats import rankdata
NP=2000
out={}
def perm_generic(f,y,seed,np_=NP):
    yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    for j in range(np_):
        yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    return dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),**mcc_conf(y,yh)), yh
# 1) observed nested yhat for A1, A1u, A2 -> bootstrap CI + prevalence threshold
obs={}
for k in KEYS:
    y=Y[k]; d=D[k]
    _,_,yh1=ao.alt_nested_cv(d,G,y,1); _,_,yhu=ur.unsup_nested_cv(d,REPS_U,y,1); _,_,yh2=ng.nested_cv(d,y,1)
    obs[k]={}
    for nm,yh in [('A1',yh1),('A1u',yhu),('A2',yh2)]:
        rng=np.random.default_rng(8302); i1=np.where(y==1)[0]; i0=np.where(y==0)[0]; bq=[]
        for b in range(2000):
            idx=np.concatenate([rng.choice(i1,len(i1)),rng.choice(i0,len(i0))]); bq.append(q2(y[idx],yh[idx]))
        obs[k][nm]=dict(q2=q2(y,yh),auc=roc_auc_score(y,yh),ci95=list(np.percentile(bq,[2.5,97.5])),mcc05=mcc_conf(y,yh)['mcc'],mcc_prev=mcc_conf(y,yh,thr=y.mean())['mcc'])
    print(k,{nm:(round(v['q2'],3),np.round(v['ci95'],2).tolist(),round(v['mcc05'],2),round(v['mcc_prev'],2)) for nm,v in obs[k].items()},flush=True)
out['nested_obs']=obs; dump('E_misc',out)
# 2) Route B 25 with 1,8-cineole NumStereoCenters corrected to 0
L25=sorted(set(core.LIT25)); assert len(L25)==25
X25=core.X_of(L25); ci=core.COMPOUNDS.index('1,8-Cineole'); j=L25.index('NumStereoCenters')
X25c=X25.copy(); print('cineole NumStereo old',X25c[ci,j]); X25c[ci,j]=0
S={'HI':'H.I.','SA':'S.A.','SPneu':'S.Pneu','SPyo':'S.Pyo','PA':'P.A.'}
out['B25_cineole0']={}
for k in KEYS:
    y=core.y_of(S[k]); r,_=perm_generic(lambda yy: loo_pls(X25c,yy,1),y,702); r0=q2(y,loo_pls(X25,y,1))
    out['B25_cineole0'][k]=dict(**r,q2_as_published_values=r0); print('B25c',k,round(r0,3),'->',round(r['q2'],3),'p',r['p_q2'],'AUC',round(r['auc'],3),'p',r['p_auc'],flush=True)
dump('E_misc',out)
# 3) Morgan2: no selection PLS-DA; nested univariate filter + PLS-DA
MB=np.load('/home/claude/rev/out/morgan2_bits.npy')
def nested_filter(Xm,y):
    yh=np.zeros(N)
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False; ytr=y[tr]; R=rankdata(Xm[tr],axis=0)
        eff,pv=ng.eff_filter(R,ytr,30); keep=np.where((eff>0.3)&(pv<=0.05))[0]
        if len(keep)==0: yh[i]=ytr.mean(); continue
        Xt=Xm[tr][:,keep]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1); ok=sd>1e-6; Xt=Xt[:,ok]; mu=mu[ok]; sd=sd[ok]
        if Xt.shape[1]==0: yh[i]=ytr.mean(); continue
        m=PLSRegression(1,scale=False).fit((Xt-mu)/sd,ytr); yh[i]=m.predict((Xm[i:i+1][:,keep][:,ok]-mu)/sd).ravel()[0]
    return yh
out['morgan2']={'n_bits':int(MB.shape[1])}
for k in KEYS:
    y=Y[k]
    r1,_=perm_generic(lambda yy: loo_pls(MB,yy,1),y,8601)
    r2,_=perm_generic(lambda yy: nested_filter(MB,yy),y,8602)
    out['morgan2'][k]=dict(no_selection=r1,nested_filter=r2)
    print('Morgan2',k,'noSel Q2 %.3f p %.4f AUC %.3f p %.4f | nestedFilter Q2 %.3f p %.4f AUC %.3f p %.4f'%(r1['q2'],r1['p_q2'],r1['auc'],r1['p_auc'],r2['q2'],r2['p_q2'],r2['auc'],r2['p_auc']),flush=True)
    dump('E_misc',out)
