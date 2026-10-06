# KORRIGIERTE SEEDS: die erste Fassung dieses Laufs hat fuer alle Routen den Seed 702
# benutzt; 702 ist aber der Seed von Route B. Die veroeffentlichten Laeufe benutzen
# 31 (A1, prun_alt.py), 47 (A1u, prun_uns.py) und 7 (A2, prun.py); die Zuordnung ist
# mit check_seeds.py gegen die ersten Nullwerte der abgelegten Dateien geprueft
# (Abweichung 0). Nur mit demselben Seed ist die Differenz der p-Werte die Wirkung der
# Korrektur und nicht zusaetzliches Monte-Carlo-Rauschen.
#
# Eigenbefund: Rueckfallvorhersage ytr.mean() im leeren-Auswahl-Zweig verraet das ausgelassene Label.
# Hier: Permutations-p-Werte mit neutraler Konstante 0.5 statt ytr.mean(), sonst identische Pipeline.
import sys,os,json,time
sys.path.insert(0,'/home/claude/rev/fixcode')      # korrigierte Kopie zuerst
sys.path.insert(1,'/home/claude/rev/QSAR_selection_leakage/code')
exec(open('/home/claude/rev/common.py').read())     # setzt sys.path erneut
for _m in ['nested_generic','altorder','unsup_rep','core']:
    sys.modules.pop(_m,None)
sys.path.insert(0,'/home/claude/rev/fixcode')
import nested_generic as ng, altorder as ao, unsup_rep as ur
for m in (ng,ao,ur):
    assert '/home/claude/rev/fixcode' in m.__file__, m.__file__
from scipy.stats import rankdata
NP=int(os.environ.get('NP','2000'))
ROUTES=os.environ.get('ROUTES','A1u,A1,A2').split(',')
SEED={'A1':31,'A1u':47,'A2':7}   # wie in den veroeffentlichten Laeufen
FOLDS=[]
for i in range(N):
    tr=np.ones(N,bool); tr[i]=False
    FOLDS.append(dict(tr=tr,i=i,Xtr=X[tr],R=rankdata(X[tr],axis=0),Xte=X[i:i+1]))
out={}
p=os.path.join('/home/claude/rev/out','I_fallback_seeds.json')
if os.path.exists(p): out=json.load(open(p))
out.setdefault('meta',dict(nperm=NP,seed={'A1':31,'A1u':47,'A2':7},change='empty-selection fallback prediction changed from mean of the training labels to the constant 0.5'))
for k in KEYS:
    y=Y[k]; D1=dict(DESC=DESC,X=X,Y=y,N=N,FOLDS=FOLDS)
    out.setdefault(k,{})
    fns={'A1':lambda yy: ao.alt_nested_cv(D1,G,yy,1)[2],
         'A1u':lambda yy: ur.unsup_nested_cv(D1,REPS_U,yy,1)[2],
         'A2':lambda yy: ng.nested_cv(D1,yy,1)[2]}
    for nm in ROUTES:
        if nm in out[k]: print('cache',k,nm,flush=True); continue
        f=fns[nm]; t=time.time()
        yh=f(y); oq=q2(y,yh); oa=roc_auc_score(y,yh); qs=[];as_=[]
        for j in range(NP):
            yp=np.random.default_rng([SEED[nm],j]).permutation(y); h=f(yp)
            qs.append(q2(yp,h)); as_.append(roc_auc_score(yp,h))
        out[k][nm]=dict(q2=oq,auc=oa,p_q2=pval(qs,oq),p_auc=pval(as_,oa),nperm=NP,seconds=round(time.time()-t,1))
        print('%s %-4s Q2 %6.3f p %.4f AUC %.3f p %.4f  %.0fs'%(k,nm,oq,out[k][nm]['p_q2'],oa,out[k][nm]['p_auc'],time.time()-t),flush=True)
        json.dump(out,open(p,'w'),indent=1,default=float)
print('fertig',ROUTES,flush=True)
