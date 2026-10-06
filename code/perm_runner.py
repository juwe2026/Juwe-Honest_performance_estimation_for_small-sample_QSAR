# Revision 03.10.2026: Permutations-p-Werte zaehlen Gleichstaende mit (Toleranz 1e-9), siehe CHANGELOG_revision.md
import sys, os, numpy as np, nested_h as nh

STORE = '/home/claude/perm_HI.npz'
NLV = 1
batch = int(sys.argv[1]) if len(sys.argv) > 1 else 400

# Beobachtet einmalig speichern
if os.path.exists(STORE):
    d = np.load(STORE)
    q2s = list(d['q2s']); aucs = list(d['aucs'])
    obs_q2 = float(d['obs_q2']); obs_auc = float(d['obs_auc'])
else:
    obs_q2, obs_auc, _, nf = nh.nested_cv(nh.Y, nlv=NLV)
    q2s, aucs = [], []
    print(f'Beobachtet: Q2={obs_q2:.4f} AUC={obs_auc:.4f}')

start = len(q2s)
for k in range(start, start + batch):
    rng = np.random.default_rng([42, k])          # reproduzierbar, unabhaengig je Index
    yv = rng.permutation(nh.Y)
    q, a, _, _ = nh.nested_cv(yv, nlv=NLV)
    q2s.append(q); aucs.append(a)

np.savez(STORE, q2s=np.array(q2s), aucs=np.array(aucs),
         obs_q2=obs_q2, obs_auc=obs_auc)
n = len(q2s)
pQ = (1 + np.sum(np.array(q2s) >= obs_q2 - 1e-9)) / (1 + n)
pA = (1 + np.sum(np.array(aucs) >= obs_auc - 1e-9)) / (1 + n)
print(f'Perms gesamt: {n} | p(Q2)={pQ:.4f} p(AUC)={pA:.4f} | '
      f'perm-mean Q2={np.mean(q2s):.3f} AUC={np.mean(aucs):.3f}')
