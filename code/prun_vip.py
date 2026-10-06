# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
# Revision 05.10.2026: ruft die rekonstruierte Funktion a2vip.nested_cv_a2_viptop auf, setzt die Aktivitaetsspalte
# ausdruecklich und schreibt archivrelativ (CHANGELOG_revision.md, Abschnitt 5).
"""Whole-pipeline permutation test for the Route A2 top-15-by-VIP arm (Additional file 10).

Usage:  python code/prun_vip.py SHORT SHEET NPERM CODE [BUDGET]

  SHORT  strain key, one of HI, SA, SPneu, SPyo, PA (also selects the activity column)
  SHEET  sheet of Additional file 6 holding the descriptor pool: Cleaned_Steps1-4
  NPERM  number of permutations (2000 in the recomputation, 1,900 to 1,982 as deposited)
  CODE   file suffix for the output, conventionally the same as SHORT
  BUDGET optional wall-clock limit in seconds. The run resumes where it stopped.

Example:  python code/prun_vip.py HI Cleaned_Steps1-4 2000 HI 100000

The script as submitted called nested_generic.nested_cv_viptop, a function missing from the
deposited module. It now calls the reconstruction a2vip.nested_cv_a2_viptop. The VIP ranking
uses two PLS components (nlv_vip=2, the default), as in the original call. Seed: permutation k
uses numpy.random.default_rng([23, k]). Checked against the archived
results/raw_permutations/perm_vip_*.npz for all five strains: the observed Q2 and AUC agree
exactly, and so do the first 20 null values of each strain, with two exceptions (P. aeruginosa,
permutations 11 and 16). In each of them one training fold reaches the fallback branch, which
now predicts 0.5 instead of the training-class mean (CHANGELOG_revision.md, section 1). With the
former fallback both values agree exactly as well.
The output goes to results/recalc/perm_vip_CODE.npz (QSAR_RECALC_DIR overrides the folder).
"""
import sys, os, time, numpy as np, nested_generic as ng, a2vip
short=sys.argv[1]; sheet=sys.argv[2]; target=int(sys.argv[3]); code=sys.argv[4]
budget=float(sys.argv[5]) if len(sys.argv)>5 else 250.0
# The pooled sheet Cleaned_Steps1-4 carries all five activity columns side by side, so the
# strain has to be named explicitly (same mapping as prun_a1.py and the README).
YCOL={'HI':'H. I. (H. influenzae)','SA':'S. A. (S. aureus)',
      'SPneu':'S. Pneu (S. pneumoniae)','SPyo':'S. pyo (S. pyogenes)',
      'PA':'P. A. (P. aeruginosa)'}
assert short in YCOL, "SHORT must be one of %s, got %r" % (sorted(YCOL), short)
STORE=os.path.join(ng._RES, f'perm_vip_{code}.npz'); t0=time.time()
D=ng.build(sheet, ycol=YCOL[short] if sheet=='Cleaned_Steps1-4' else None)
if os.path.exists(STORE):
    d=np.load(STORE); q2s=list(d['q2s']); aucs=list(d['aucs']); oq=float(d['obs_q2']); oa=float(d['obs_auc'])
else:
    oq,oa,_=a2vip.nested_cv_a2_viptop(D,D['Y'],1,15); q2s=[]; aucs=[]
def save(): np.savez(STORE,q2s=np.array(q2s),aucs=np.array(aucs),obs_q2=oq,obs_auc=oa)
while len(q2s)<target and (time.time()-t0)<budget:
    k=len(q2s); yv=np.random.default_rng([23,k]).permutation(D['Y'])
    q,a,_=a2vip.nested_cv_a2_viptop(D,yv,1,15); q2s.append(q); aucs.append(a)
    if len(q2s)%50==0: save()
save(); n=len(q2s)
pQ=(1+np.sum(np.array(q2s)>=oq - 1e-9))/(1+n); pA=(1+np.sum(np.array(aucs)>=oa - 1e-9))/(1+n)
print(f'{short} VIP-top15: {n}/{target} | obs Q2={oq:.3f} AUC={oa:.3f} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | {time.time()-t0:.0f}s')
