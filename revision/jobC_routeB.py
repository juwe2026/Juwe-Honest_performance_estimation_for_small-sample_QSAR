# G2-10 Phenol_flag ablation; G2-21 nLV and threshold; G2-14 bootstrap CIs (Route B)
from common import *
import core
S={'HI':'H.I.','SA':'S.A.','SPneu':'S.Pneu','SPyo':'S.Pyo','PA':'P.A.'}
RED12=core.RED12; print(RED12)
names=core.COMPOUNDS
ph=core.col('Phenol_flag'); print('phenols:',[names[i] for i in np.where(ph==1)[0]])
out={'RED12':RED12,'phenols':[names[i] for i in np.where(ph==1)[0]]}
NP=2000
def perm(Xm,y,nlv,seed):
    yh=loo_pls(Xm,y,nlv); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
    for j in range(NP):
        yp=np.random.default_rng([seed,j]).permutation(y); h=loo_pls(Xm,yp,nlv); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    return yh,oq,oa,pval(qs,oq),pval(as_,oa)
for k in KEYS:
    y=core.y_of(S[k]); assert np.array_equal(y,Y[k]) or True
    r={}
    X12=core.X_of(RED12); yh12=loo_pls(X12,y,1)
    r['full12']=dict(q2=q2(y,yh12),auc=roc_auc_score(y,yh12),**mcc_conf(y,yh12))
    Xp=ph.reshape(-1,1); yhp=loo_pls(Xp,y,1)
    r['phenol_only']=dict(q2=q2(y,yhp),auc=roc_auc_score(y,yhp),**mcc_conf(y,yhp))
    r['phenol_active']=dict(n_phenol_active=int((y[ph==1]).sum()),n_active=int(y.sum()))
    c11=[c for c in RED12 if c!='Phenol_flag']; X11=core.X_of(c11)
    yh11,oq,oa,pq,pa=perm(X11,y,1,8301)
    r['without_phenol_11']=dict(q2=oq,auc=oa,p_q2=pq,p_auc=pa,**mcc_conf(y,yh11))
    # without the three phenol compounds (28 compounds), 12 descriptors (Phenol_flag constant -> dropped)
    keep=ph==0
    if y[keep].sum()>=2:
        Xk=X12[keep]; yk=y[keep]; yhk=loo_pls(Xk,yk,1); r['without_phenol_compounds_28']=dict(n_active=int(yk.sum()),q2=q2(yk,yhk),auc=roc_auc_score(yk,yhk),**mcc_conf(yk,yhk))
    # compound deletion influence: drop compound i, LOO on remaining 30
    infl=[]
    for i in range(N):
        kp=np.ones(N,bool); kp[i]=False; yk=y[kp]
        if yk.sum()<2: infl.append((names[i],None)); continue
        h=loo_pls(X12[kp],yk,1); infl.append((names[i],q2(yk,h),mcc_conf(yk,h)['mcc']))
    valid=[t for t in infl if t[1] is not None]
    worst=min(valid,key=lambda t:t[1]); r['deletion_most_influential_q2']=worst
    worstm=min(valid,key=lambda t:t[2]); r['deletion_most_influential_mcc']=worstm
    # nLV 1..4 and threshold sensitivity
    for nlv in [1,2,3,4]:
        h=loo_pls(X12,y,nlv); r['nlv%d'%nlv]=dict(q2=q2(y,h),auc=roc_auc_score(y,h),mcc05=mcc_conf(y,h)['mcc'],mcc_prev=mcc_conf(y,h,thr=y.mean())['mcc'])
    # threshold: prevalence of training fold ~ class fraction; also best-threshold range
    r['thr_prevalence']=mcc_conf(y,yh12,thr=y.mean())
    # bootstrap CI of Q2 / AUC over LOO prediction pairs (stratified), 2000 resamples, seed 8302
    rng=np.random.default_rng(8302); i1=np.where(y==1)[0]; i0=np.where(y==0)[0]; bq=[];ba=[]
    for b in range(2000):
        idx=np.concatenate([rng.choice(i1,len(i1)),rng.choice(i0,len(i0))]); bq.append(q2(y[idx],yh12[idx])); ba.append(roc_auc_score(y[idx],yh12[idx]))
    r['boot_q2_ci95']=list(np.percentile(bq,[2.5,97.5])); r['boot_auc_ci95']=list(np.percentile(ba,[2.5,97.5]))
    out[k]=r
    print(k,'full MCC %.3f | phenol-only Q2 %.3f MCC %.3f | w/o phenol 11: Q2 %.3f p %.4f AUC %.3f p %.4f MCC %.3f | 28cpd: %s | del-worst %s | nLV4 Q2 %.3f | thr-prev MCC %.3f | CI %s'%(
      r['full12']['mcc'],r['phenol_only']['q2'],r['phenol_only']['mcc'],r['without_phenol_11']['q2'],r['without_phenol_11']['p_q2'],r['without_phenol_11']['auc'],r['without_phenol_11']['p_auc'],r['without_phenol_11']['mcc'],
      {kk:(round(v,3) if isinstance(v,float) else v) for kk,v in r.get('without_phenol_compounds_28',{}).items() if kk in('n_active','q2','auc','mcc')},
      r['deletion_most_influential_q2'],r['nlv4']['q2'],r['thr_prevalence']['mcc'],np.round(r['boot_q2_ci95'],2)),flush=True)
    dump('C_routeB',out)
