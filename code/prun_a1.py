import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 4 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import sys, os, time, numpy as np, nested_generic as ng, altorder as ao
"""Whole-pipeline permutation test for the Route A1 top-15 models.

Usage:  python code/prun_a1.py MODE SHORT SHEET NPERM CODE [BUDGET]

  MODE   e = top 15 by effect size (Fig. 5, Table 3), v = top 15 by VIP (Table 4)
  SHORT  strain key, one of HI, SA, SPneu, SPyo, PA (also selects the activity column)
  SHEET  sheet of Additional file 4 holding the descriptor pool: Cleaned_Steps1-4
  NPERM  number of permutations (2000 in the published runs)
  CODE   file suffix for the output, conventionally the same as SHORT
  BUDGET optional wall-clock limit in seconds; the run resumes where it stopped

Example:  python code/prun_a1.py e HI Cleaned_Steps1-4 2000 HI 100000
"""
mode=sys.argv[1]; short=sys.argv[2]; sheet=sys.argv[3]; target=int(sys.argv[4]); code=sys.argv[5]
budget=float(sys.argv[6]) if len(sys.argv)>6 else 250.0
# The pooled sheet Cleaned_Steps1-4 carries all five activity columns side by side, so
# the strain has to be named explicitly. Falling back on the third column, as the
# strain-specific CS1-5_* sheets allow, would silently run every strain against
# H. influenzae's labels (round 6 report, 9.5).
YCOL={'HI':'H. I. (H. influenzae)','SA':'S. A. (S. aureus)',
      'SPneu':'S. Pneu (S. pneumoniae)','SPyo':'S. pyo (S. pyogenes)',
      'PA':'P. A. (P. aeruginosa)'}
assert short in YCOL, "SHORT must be one of %s, got %r" % (sorted(YCOL), short)
STORE=_os.path.join(_RES, f'perm_a1{mode}_{code}.npz'); t0=time.time()
D=ng.build(sheet, ycol=YCOL[short] if sheet=='Cleaned_Steps1-4' else None); groups=ao.get_groups(D['X'])
fn = (lambda yv: ao.alt_nested_cv(D,groups,yv,1,top_k=15)) if mode=='e' else (lambda yv: ao.alt_viptop(D,groups,yv,1,15))
if os.path.exists(STORE):
    d=np.load(STORE); q2s=list(d['q2s']); aucs=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
else:
    oq,oa,_=fn(D['Y']); q2s=[]; aucs=[]
def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
while len(q2s)<target and (time.time()-t0)<budget:
    k=len(q2s); yv=np.random.default_rng([59 if mode=='e' else 73,k]).permutation(D['Y'])
    q,a,_=fn(yv); q2s.append(q); aucs.append(a)
    if len(q2s)%200==0: save()
save(); n=len(q2s)
pQ=(1+np.sum(np.array(q2s)>=oq))/(1+n); pA=(1+np.sum(np.array(aucs)>=oa))/(1+n)
print(f'{short} A1-{"eff" if mode=="e" else "vip"}15: {n}/{target} | Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | {time.time()-t0:.0f}s')
