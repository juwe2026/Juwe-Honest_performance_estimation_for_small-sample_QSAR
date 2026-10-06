# KORRIGIERTE FASSUNG des Morgan-Arms mit genestetem Filter (jobE_misc.py, Teil 3).
# Die eigene Filterfunktion setzte bei leerer Auswahl ytr.mean(); bei den Morgan-Bits
# passiert fuer H. influenzae und S. pyogenes in den meisten Folds kein Bit den Filter,
# der Zweig greift also schon bei den echten Labels und verfaelscht den Punktschaetzer
# (AUC 0,000 = perfekt umgekehrte Rangfolge). Jetzt die neutrale Konstante 0.5,
# gleicher Seed 8602, 2.000 Permutationen. Der Arm ohne Auswahl hat keinen solchen
# Zweig und bleibt unveraendert.
import sys, os
os.environ.setdefault('QSAR_ROUTEA_XLSX','/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx')
sys.path.insert(0,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())
import core
for _m in ['nested_generic','altorder','unsup_rep']:
    sys.modules.pop(_m,None)
sys.path.insert(0,'/home/claude/rev/fixcode')
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur): assert '/home/claude/rev/fixcode' in m.__file__, m.__file__
from scipy.stats import rankdata
import json
NP=2000
MB=np.load('/home/claude/rev/out/morgan2_bits.npy')
FB=[0]
def nested_filter(Xm,y):
    yh=np.zeros(N)
    for i in range(N):
        tr=np.ones(N,bool); tr[i]=False; ytr=y[tr]; R=rankdata(Xm[tr],axis=0)
        eff,pv=ng.eff_filter(R,ytr,30); keep=np.where((eff>0.3)&(pv<=0.05))[0]
        if len(keep)==0: yh[i]=0.5; FB[0]+=1; continue
        Xt=Xm[tr][:,keep]; mu=Xt.mean(0); sd=Xt.std(0,ddof=1); ok=sd>1e-6; Xt=Xt[:,ok]; mu=mu[ok]; sd=sd[ok]
        if Xt.shape[1]==0: yh[i]=0.5; FB[0]+=1; continue
        m=PLSRegression(1,scale=False).fit((Xt-mu)/sd,ytr); yh[i]=m.predict((Xm[i:i+1][:,keep][:,ok]-mu)/sd).ravel()[0]
    return yh
E=json.load(open('/home/claude/rev/out/E_misc.json'))
out={'n_bits':int(MB.shape[1])}
for k in KEYS:
    y=Y[k]; FB[0]=0
    yh=nested_filter(MB,y); obs_fb=FB[0]; oq=q2(y,yh)
    try: oa=roc_auc_score(y,yh)
    except Exception: oa=float('nan')
    qs=[];as_=[]
    for j in range(NP):
        yp=np.random.default_rng([8602,j]).permutation(y); h=nested_filter(MB,yp)
        qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
    old=E['morgan2'][k]['nested_filter']
    out[k]=dict(no_selection=E['morgan2'][k]['no_selection'],
                nested_filter=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),**mcc_conf(y,yh),
                                   empty_folds_observed=obs_fb,q2_old=old['q2'],auc_old=old['auc'],
                                   p_q2_old=old['p_q2'],p_auc_old=old['p_auc']))
    r=out[k]['nested_filter']
    print('%-6s leere Folds %2d/31 | Q2 %.3f (alt %.3f) p %.4f | AUC %.3f (alt %.3f) p %.4f'%(
        k,obs_fb,oq,old['q2'],r['p_q2'],oa,old['auc'],r['p_auc']),flush=True)
    json.dump(out,open('/home/claude/rev/out/morgan2_fixed.json','w'),indent=1,default=float)
print('fertig',flush=True)
