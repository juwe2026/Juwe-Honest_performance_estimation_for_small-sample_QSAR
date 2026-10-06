# G2-14: Route A1 / A1u without the conformer-dependent 3D descriptor blocks (Chemopy_3D, PaDEL_3D)
# KORRIGIERTE FASSUNG: korrigierter Rueckfallzweig (fixcode), gleiche Seeds 8701/8702.
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
keep=np.array([not (d.startswith('Chemopy_3D::') or d.startswith('PaDEL_3D::')) for d in DESC])
X2=X[:,keep]; print('2D pool',X2.shape)
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False; FOLDS.append(dict(tr=tr,i=i,Xtr=X2[tr],R=rankdata(X2[tr],axis=0),Xte=X2[i:i+1]))
G2,_=ao.build_groups(X2); R2=ur.centrality_reps(X2,G2); print('groups',len(G2))
NP=2000; out={'n_desc':int(X2.shape[1]),'n_groups':len(G2)}
for k in KEYS:
    y=Y[k]; Dk=dict(DESC=None,X=X2,Y=y,N=N,FOLDS=FOLDS); out[k]={}
    for nm,f,seed in [('A1',lambda yy: ao.alt_nested_cv(Dk,G2,yy,1)[2],8701),('A1u',lambda yy: ur.unsup_nested_cv(Dk,R2,yy,1)[2],8702)]:
        t=time.time(); yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
        for j in range(NP):
            yp=np.random.default_rng([seed,j]).permutation(y); h=f(yp); qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
        out[k][nm]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa))
        print(k,nm,'Q2 %.3f p %.4f AUC %.3f p %.4f %.0fs'%(oq,pval(qs,oq),oa,pval(as_,oa),time.time()-t),flush=True)
        dump('F_no3d_fixed',out)
