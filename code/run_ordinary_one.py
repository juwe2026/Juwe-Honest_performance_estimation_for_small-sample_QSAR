import sys, pickle, os, numpy as np, nested_generic as ng, altorder as ao, unsup_rep as ur, ordinary_perm as op
from scipy.stats import rankdata
short = sys.argv[1]; sheet = sys.argv[2]
STORE = '/home/claude/ordinary_perm_results.pkl'
results = pickle.load(open(STORE,'rb')) if os.path.exists(STORE) else {}
if short in results:
    print(f"{short}: bereits vorhanden, ueberspringe."); sys.exit(0)
D0 = ng.build('CS1-4_H.I.'); groups = ao.get_groups(D0['X']); reps_u = ur.centrality_reps(D0['X'], groups)
D = ng.build(sheet); y = D['Y']; N = D['N']; R = rankdata(D['X'], axis=0)
cols1, eff1, pv1, _ = ao.alt_select(groups, R, y, N)
r1 = op.ordinary_perm_test(D['X'], y, cols1, 1, n_perm=1000, seed=201)
ecomb, pcomb = ng.eff_filter(R, y, N); keep = np.where((ecomb>0.3)&(pcomb<=0.05))[0]
rl = ng.group_reps(D['X'][:,keep], R[:,keep], ecomb[keep], N); cols2 = keep[rl]
r2 = op.ordinary_perm_test(D['X'], y, cols2, 1, n_perm=1000, seed=202)
colsv = ng._viptop(D['X'], y, cols1, 15) if len(cols1) > 15 else cols1
rv = op.ordinary_perm_test(D['X'], y, colsv, 1, n_perm=1000, seed=203)
colsu, effu, pvu = ur.unsup_select(reps_u, R, y, N)
ru = op.ordinary_perm_test(D['X'], y, colsu, 1, n_perm=1000, seed=204)
results[short] = dict(A1=r1[:4], A2=r2[:4], VIP=rv[:4], A1u=ru[:4])
pickle.dump(results, open(STORE, 'wb'))
print(f"{short} A1 q={r1[0]:+.3f}(p={r1[2]:.3f}) a={r1[1]:.3f}(p={r1[3]:.3f}) | "
      f"A2 q={r2[0]:+.3f}(p={r2[2]:.3f}) a={r2[1]:.3f}(p={r2[3]:.3f}) | "
      f"VIP q={rv[0]:+.3f}(p={rv[2]:.3f}) a={rv[1]:.3f}(p={rv[3]:.3f}) | "
      f"A1u q={ru[0]:+.3f}(p={ru[2]:.3f}) a={ru[1]:.3f}(p={ru[3]:.3f})")
