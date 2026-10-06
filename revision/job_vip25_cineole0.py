# -*- coding: utf-8 -*-
"""G2-20: VIP (zwei Komponenten) der 25 Literaturdeskriptoren, apparent und Mittel ueber die LOO-Folds, mit Cineol = 2
(Reproduktion der eingereichten Werte) und = 0. Ausgabe out/vip25_cineole0.json (Figs. S2c/S3b, VIP-Tabellen)."""
import sys; sys.path.insert(0,'/home/claude/rev')
from common import *
import core, numpy as np, pandas as pd
from sklearn.cross_decomposition import PLSRegression
L25=sorted(set(core.LIT25)); X=core.X_of(L25); ci=core.COMPOUNDS.index('1,8-Cineole'); j=L25.index('NumStereoCenters')
X0=X.copy(); X0[ci,j]=0
def vipfit(Xr,y,nlv=2):
    Xs=core.autoscale_all(Xr); m=PLSRegression(n_components=nlv,scale=False).fit(Xs,y); return core.vip(m,Xs)
def nested(Xr,y,nlv=2):
    V=[]
    for i in range(len(y)):
        tr=np.ones(len(y),bool); tr[i]=False; V.append(vipfit(Xr[tr],y[tr],nlv))
    V=np.array(V); return V.mean(0), (V>1).sum(0)
res={}
for s in core.STRAINS:
    y=core.y_of(s); W=pd.read_excel('/home/claude/rev/QSAR_selection_leakage/data/20260710_VIP_EffectSize_Tabellen.xlsx','VIP_Lit_'+s).set_index('Descriptor')
    a2=vipfit(X,y); n2,f2=nested(X,y); a0=vipfit(X0,y); n0,f0=nested(X0,y)
    da=max(abs(W.loc[d,'TwoStage_VIP_apparent']-a2[k]) for k,d in enumerate(L25)); dn=max(abs(W.loc[d,'TwoStage_VIP_nested_mean']-n2[k]) for k,d in enumerate(L25))
    fq=[str(W.loc[d,'TwoStage_selfreq_VIPgt1'])=='%d/31'%f2[k] for k,d in enumerate(L25)]
    print(s,'reproduce apparent %.4f nested %.4f selfreq all %s'%(da,dn,all(fq)), '| change max excl NumStereo app %.3f nest %.3f'%(np.delete(abs(a0-a2),j).max(),np.delete(abs(n0-n2),j).max()),
          '| NumStereo app %.3f->%.3f nest %.3f->%.3f freq %d->%d'%(a2[j],a0[j],n2[j],n0[j],f2[j],f0[j]))
    res[s]=dict(app=a0.tolist(),nest=n0.tolist(),freq=f0.tolist(),app2=a2.tolist(),nest2=n2.tolist(),freq2=f2.tolist())
import json; json.dump(dict(L25=L25,res=res),open('/home/claude/rev/out/vip25_cineole0.json','w'))
