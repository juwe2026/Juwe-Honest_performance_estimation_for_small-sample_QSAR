# Effective number of independent endpoints (Li & Ji 2005, Nyholt 2004 / Cheverud 2001) for AF11 and Table S20 of Revision 1
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)
exec(open((_REV + '/common.py')).read())
import numpy as np, json
import core; names=list(core.COMPOUNDS)
RAW=json.load(open((_OUT + '/raw_mic_31.json')))
# Check: classes = mean rule from the raw values, order of the compounds
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
R=json.load(open((_REV + '/lone/rows.json')))
THR={'A':{'EC':8,'SA':8,'PA':16,'AB':8},'C':{'EC':16,'SA':16,'PA':16,'AB':16}}
def mv(x):
    x=str(x).strip(); return None if x.startswith('>') else float(x)
for r in THR:
    L2=np.array([[int(mv(x['MIC_'+k]) is not None and mv(x['MIC_'+k])<=THR[r][k]) for x in R] for k in THR[r]],float)
    C2,e2,lj2,ny2=meff(L2)
    out['lone_rule_'+r]=dict(strains=list(THR[r]),phi=C2.tolist(),eigen=e2.tolist(),li_ji=lj2,nyholt=ny2)
json.dump(out,open((_OUT + '/endpoint_corr.json'),'w'),indent=1)
print('mono eigenvalues',np.round(ev,3),'LiJi %.2f Nyholt %.2f'%(lj,ny))
print('phi min/max',np.round(C[np.triu_indices(5,1)].min(),3),np.round(C[np.triu_indices(5,1)].max(),3),'both',min(both.values()),max(both.values()))
for r in THR: print(r,'LiJi %.2f Nyholt %.2f'%(out['lone_rule_'+r]['li_ji'],out['lone_rule_'+r]['nyholt']))
