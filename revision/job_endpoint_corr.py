# Effektive Zahl unabhaengiger Endpunkte (Li & Ji 2005; Nyholt 2004 / Cheverud 2001) fuer AF11 und Table S20
exec(open('/home/claude/rev/common.py').read())
import numpy as np, json
import core; names=list(core.COMPOUNDS)
RAW=json.load(open('/home/claude/rev/out/raw_mic_31.json'))
# Pruefung: Klassen = Mittelwert-Regel aus den Rohwerten, Reihenfolge der Verbindungen
assert names is not None and len(names)==31
for k in KEYS:
    for nm,y in zip(names,Y[k]):
        m=RAW[nm][k+'_broth']['mean']
        assert int(m<2048)==int(y),(k,nm,m,y)
def meff(L):
    C=np.corrcoef(L); ev=np.sort(np.linalg.eigvalsh(C))[::-1]; m=len(ev)
    lj=float(sum((e>=1)+(e-np.floor(e)) for e in np.abs(ev)))
    ny=float(1+(m-1)*(1-np.var(ev,ddof=1)/m))
    return C,ev,lj,ny
L=np.array([np.asarray(Y[k],float) for k in KEYS])
C,ev,lj,ny=meff(L)
both={'%s|%s'%(a,b):int(((L[i]==1)&(L[j]==1)).sum()) for i,a in enumerate(KEYS) for j,b in enumerate(KEYS) if i<j}
out=dict(mono=dict(strains=KEYS,compounds=names,Y={k:[int(v) for v in Y[k]] for k in KEYS},phi=C.tolist(),eigen=ev.tolist(),
                   li_ji=lj,nyholt=ny,both_active=both,n_active={k:int(np.sum(Y[k])) for k in KEYS}))
R=json.load(open('/home/claude/rev/lone/rows.json'))
THR={'A':{'EC':8,'SA':8,'PA':16,'AB':8},'C':{'EC':16,'SA':16,'PA':16,'AB':16}}
def mv(x):
    x=str(x).strip(); return None if x.startswith('>') else float(x)
for r in THR:
    L2=np.array([[int(mv(x['MIC_'+k]) is not None and mv(x['MIC_'+k])<=THR[r][k]) for x in R] for k in THR[r]],float)
    C2,e2,lj2,ny2=meff(L2)
    out['lone_rule_'+r]=dict(strains=list(THR[r]),phi=C2.tolist(),eigen=e2.tolist(),li_ji=lj2,nyholt=ny2)
json.dump(out,open('/home/claude/rev/out/endpoint_corr.json','w'),indent=1)
print('mono eigen',np.round(ev,3),'LiJi %.2f Nyholt %.2f'%(lj,ny))
print('phi min/max',np.round(C[np.triu_indices(5,1)].min(),3),np.round(C[np.triu_indices(5,1)].max(),3),'both',min(both.values()),max(both.values()))
for r in THR: print(r,'LiJi %.2f Nyholt %.2f'%(out['lone_rule_'+r]['li_ji'],out['lone_rule_'+r]['nyholt']))
