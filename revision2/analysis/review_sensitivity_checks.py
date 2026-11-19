# -*- coding: utf-8 -*-
"""Three sensitivity checks on observed values (no permutation except for the leave-group-out check).

  sd_rule      Routes A1, A1u and A2 with and without the descriptors whose standard deviation over all compounds
               is below 1e-6 (monoterpenoids: two columns, peptides: one column). The arms without filter discard
               such columns in every fold, the three routes do not.
  rf_seeds     out-of-bag Q2 and AUC of the random forest of the peptides (rule A) for the seeds 0 to 19.
               Table 2 uses seed 0.
  routeB_lgo   Route B (12 descriptors, one latent variable): leave-one-out against leave-group-out, where
               compounds with an identical descriptor vector are left out together, 2,000 permutations (seed 9401).

Usage:  python review_sensitivity_checks.py <archive_root> <out.json>
Needs Additional file 3 in data/ (or QSAR_ROUTEA_XLSX) and an empty QSAR_RECALC_DIR. Running time a few minutes.
"""
import sys, os, json, collections
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestClassifier
R, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(R, 'code'))
import nested_generic as ng, altorder as ao, unsup_rep as ur
def q2(y, yh): return float(1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2))
def auc(y, yh): return float(roc_auc_score(y, yh))
res = {}

def routes(X, DESC, Y):
    N = X.shape[0]; FOLDS = []
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        FOLDS.append(dict(tr=tr, i=i, Xtr=X[tr], R=rankdata(X[tr], axis=0), Xte=X[i:i + 1]))
    G, _ = ao.build_groups(X); REPS = np.asarray(ur.centrality_reps(X, G)); o = {'groups': len(G)}
    for k, y in Y.items():
        D = dict(DESC=DESC, X=X, Y=y, N=N, FOLDS=FOLDS)
        o[k] = {nm: [round(q2(y, h), 4), round(auc(y, h), 4)] for nm, h in (
            ('A1', ao.alt_nested_cv(D, G, y, 1)[2]), ('A1u', ur.unsup_nested_cv(D, REPS, y, 1)[2]), ('A2', ng.nested_cv(D, y, 1)[2]))}
    return o
def sd_rule(X, DESC, Y):
    sd = X.std(axis=0); low = np.where(sd < 1e-6)[0]; keep = np.setdiff1d(np.arange(X.shape[1]), low)
    fold_sd = max(float(np.delete(X[:, j], i).std(ddof=1)) for j in low for i in range(X.shape[0]))
    return dict(columns=[DESC[j] for j in low], largest_training_sd=fold_sd, all_columns=routes(X, DESC, Y),
                without_these_columns=routes(X[:, keep], [DESC[j] for j in keep], Y))
# monoterpenoids
D0 = ng.build('Cleaned_Steps1-4', ycol='S. A. (S. aureus)')
YC = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPn': 'S. Pneu (S. pneumoniae)', 'SPy': 'S. pyo (S. pyogenes)'}
YM = {k: (ng._clean[c].astype(str).str.strip() == 'active').astype(float).values for k, c in YC.items()}
res['sd_rule'] = {'monoterpenoids': sd_rule(np.asarray(D0['X'], float), list(D0['DESC']), YM)}
# peptides
rows = json.load(open(os.path.join(R, 'revision', 'lone', 'rows.json'))); P = pd.read_csv(os.path.join(R, 'revision', 'lone', 'desc_pool.csv'))
THR = {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}
YP = {k: np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows]) for k in THR}
DP = [c for c in P.columns if c != 'No']; XP = P[DP].values.astype(float)
res['sd_rule']['peptides'] = sd_rule(XP, DP, YP)
for ds, v in res['sd_rule'].items():
    print(ds, v['columns'], 'largest training SD %.2g' % v['largest_training_sd'], 'groups', v['all_columns']['groups'], '->', v['without_these_columns']['groups'])
    for k in v['all_columns']:
        if k != 'groups': print('   %-4s all %s | without %s' % (k, v['all_columns'][k], v['without_these_columns'][k]), flush=True)
# random forest seeds
rf = {}
for k, y in YP.items():
    qs, as_ = [], []
    for seed in range(20):
        m = RandomForestClassifier(300, min_samples_leaf=1, max_features='sqrt', oob_score=True, n_jobs=1, random_state=seed).fit(XP, y)
        p = m.oob_decision_function_[:, 1]; qs.append(round(q2(y, p), 4)); as_.append(round(auc(y, p), 4))
    rf[k] = dict(q2=qs, auc=as_, q2_seed0=qs[0], q2_min=min(qs), q2_median=float(np.median(qs)), q2_max=max(qs), auc_seed0=as_[0], auc_min=min(as_), auc_median=float(np.median(as_)), auc_max=max(as_))
    print('RF %-3s Q2 seed0 %.3f range %.3f to %.3f median %.3f | AUC seed0 %.3f range %.3f to %.3f' % (k, qs[0], min(qs), max(qs), np.median(qs), as_[0], min(as_), max(as_)), flush=True)
res['rf_seeds'] = rf
# Route B, leave-group-out
RB = pd.read_excel(os.path.join(R, 'data', 'RouteB_descriptor_data.xlsx'), sheet_name='Descriptors')
RED12 = ['TPSA', 'RingCount', 'FractionCSP3', 'NumAromaticRings', 'CeqC_DoubleBonds', 'Aldehyde_flag', 'Ketone_flag', 'CarboxAcid_flag',
         'Phenol_flag', 'Ether_flag', 'EsterLacton_flag', 'Alcohol_flag']
XB = RB[RED12].to_numpy(float); NB = len(XB); gid = {}; grp = np.array([gid.setdefault(XB[i].tobytes(), len(gid)) for i in range(NB)])
def cv_pls(y, groups):
    yh = np.zeros(NB)
    for g in np.unique(groups):
        te = groups == g; tr = ~te; Xt = XB[tr]; mu = Xt.mean(0); sd = Xt.std(0, ddof=1); ok = sd > 1e-6
        m = PLSRegression(n_components=1, scale=False).fit((Xt[:, ok] - mu[ok]) / sd[ok], y[tr])
        yh[te] = m.predict((XB[te][:, ok] - mu[ok]) / sd[ok]).ravel()
    return yh
def pv(null, obs): return float((1 + np.sum(np.asarray(null) >= obs - 1e-9)) / (1 + len(null)))
lgo = {'groups': int(len(gid)), 'compounds': NB}
for c in RB.columns[5:10]:
    y = (RB[c].astype(str).str.strip() == 'active').to_numpy(float); k = c.split('(')[0].strip()
    if y.sum() < 3: continue
    o = {}
    for nm, gr in (('leave_one_out', np.arange(NB)), ('leave_group_out', grp)):
        h = cv_pls(y, gr); oq, oa = q2(y, h), auc(y, h); qs, as_ = [], []
        for j in range(2000):
            yp = np.random.default_rng([9401, j]).permutation(y); hp = cv_pls(yp, gr); qs.append(q2(yp, hp)); as_.append(auc(yp, hp))
        o[nm] = dict(q2=round(oq, 4), auc=round(oa, 4), p_q2=round(pv(qs, oq), 4), p_auc=round(pv(as_, oa), 4))
    lgo[k] = o; print('Route B %-8s LOO %s | LGO %s' % (k, o['leave_one_out'], o['leave_group_out']), flush=True)
res['routeB_lgo'] = lgo
json.dump(res, open(OUT, 'w'), indent=1)
