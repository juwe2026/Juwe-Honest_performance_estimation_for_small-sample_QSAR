# -*- coding: utf-8 -*-
"""Descriptive checks of the deposited data: missing values, replicate support of the activity classes,
identical descriptor rows, censored values and molecular size.

Usage:  python data_description_checks.py <archive_root> <Additional_file_2.xlsx> <Additional_file_3.xlsx> <out.json>
Reads only deposited files, fits no model and writes nothing into the archive. Running time a few seconds.

  missing_monoterpenoids   cells of the cleaned Route A matrix (Additional file 3) that are missing in the raw sheet
  replicates               mean rule and median rule of the 31 monoterpenoids from the replicate MIC values (Additional file 2)
  routeB_identical_rows    compounds sharing their vector of the 12 Route B descriptors
  benchmark                right-censored values and identical descriptor rows in DS01 to DS16 of the main analysis
  peptides                 class sizes under rules A and C, the six truncated peptides, AUC of the molecular weight alone,
                           columns with a standard deviation below 1e-6, cells equal to the column median for peptides 30 to 32
"""
import sys, os, json, glob, collections
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from sklearn.metrics import roc_auc_score
R, AF2, AF3, OUT = sys.argv[1:5]; res = {}

# ---- monoterpenoids: missing values
raw = pd.read_excel(AF3, sheet_name='31MT_5Bac_AllDes'); cl = pd.read_excel(AF3, sheet_name='Cleaned_Steps1-4')
meta = list(raw.columns[:7]); D = [c for c in cl.columns if c not in meta]; names = list(raw['Monoterpenoid'])
assert names == list(cl['Monoterpenoid'])
Rw = raw[D].apply(pd.to_numeric, errors='coerce'); miss = Rw.isna(); C = cl[D].to_numpy(float)
for j, c in enumerate(D):                       # cleaned value = raw value, or the median of the available values
    m = miss[c].to_numpy(); v = Rw[c].to_numpy(float)
    assert np.allclose(C[~m, j], v[~m], rtol=1e-12, atol=0) and (not m.any() or np.allclose(C[m, j], np.nanmedian(v), rtol=1e-12, atol=0))
cols = miss.sum()[miss.sum() > 0]; comp = miss.sum(axis=1).to_numpy() > 0
rawlib = [c for c in raw.columns if str(c).startswith('PaDEL_3D::')]
fis = {}
for y in meta[2:]:
    a = (raw[y].astype(str).str.strip() == 'active').to_numpy()
    tab = [[int((comp & a).sum()), int((comp & ~a).sum())], [int((~comp & a).sum()), int((~comp & ~a).sum())]]
    fis[y] = dict(table=tab, fisher_p=float(fisher_exact(tab)[1]))
sd = cl[D].std(ddof=0)
res['missing_monoterpenoids'] = dict(
    descriptors=len(D), cells=int(Rw.size), missing_cells=int(miss.sum().sum()), percent=round(100 * miss.sum().sum() / Rw.size, 2),
    imputed_columns=len(cols), imputed_columns_by_library=dict(collections.Counter(c.split('::')[0] for c in cols.index)),
    retained_by_library=dict(collections.Counter(c.split('::')[0] for c in D)),
    compounds=[n for n, k in zip(names, comp) if k], missing_per_column=sorted(set(int(v) for v in cols.values)),
    raw_PaDEL_3D_columns=len(rawlib), raw_PaDEL_3D_all_missing_for_these=bool(raw.loc[comp, rawlib].apply(pd.to_numeric, errors='coerce').isna().all().all()),
    raw_PaDEL_3D_missing_for_the_others=int(raw.loc[~comp, rawlib].apply(pd.to_numeric, errors='coerce').isna().sum().sum()),
    missing_block_versus_activity=fis, columns_sd_below_1e_6={c: float(sd[c]) for c in sd.index[sd < 1e-6]})

# ---- monoterpenoids: replicates
M = pd.read_excel(AF2, sheet_name='MIC_replicates_31'); V = pd.read_excel(AF2, sheet_name='Activity_vectors_31').set_index('Compound')
L = M[M['Phase'].str.startswith('liquid')]; g = L.groupby(['Compound', 'Strain'])['MIC_ug_mL']
n, mean, med, inh = g.count(), g.mean(), g.median(), g.apply(lambda s: int((s < 2048).sum()))
assert ((mean < 2048) == (inh > 0)).all()
def meff(Y):
    ev = np.clip(np.linalg.eigvalsh(np.corrcoef(Y, rowvar=False)), 0, None); k = len(ev)
    return dict(li_ji=round(float(sum((e >= 1) + (e - np.floor(e)) for e in ev)), 2), nyholt=round(float(1 + (k - 1) * (1 - ev.var(ddof=1) / k)), 2))
rep = {}; Ymed = {}
for s in V.columns:
    act = (mean.xs(s, level='Strain') < 2048).astype(int).reindex(V.index); assert (act == V[s]).all()
    md = (med.xs(s, level='Strain').reindex(V.index) < 2048); Ymed[s] = md.to_numpy(float)
    rep[s] = dict(active_mean_rule=int(act.sum()), active_median_rule=sorted(V.index[md]),
                  active_on_a_minority_of_replicates={c: '%d of %d' % (inh[(c, s)], n[(c, s)]) for c in V.index[act == 1] if 2 * inh[(c, s)] < n[(c, s)]})
ev4 = [s for s in V.columns if s != 'P. aeruginosa']
res['replicates'] = dict(liquid_phase_values=int(len(L)), replicates_per_compound_and_strain={int(k): int(v) for k, v in collections.Counter(n.values).items()},
    mean_below_2048_equals_at_least_one_replicate_below_2048=True, strains=rep,
    median_rule_S_pyogenes_equals_S_aureus=rep['S. pyogenes']['active_median_rule'] == rep['S. aureus']['active_median_rule'],
    effective_endpoints_four_evaluable=dict(mean_rule=meff(V[ev4].to_numpy(float)), median_rule=meff(np.column_stack([Ymed[s] for s in ev4]))))

# ---- Route B: identical descriptor rows
RB = pd.read_excel(os.path.join(R, 'data', 'RouteB_descriptor_data.xlsx'), sheet_name='Descriptors')
RED12 = ['TPSA', 'RingCount', 'FractionCSP3', 'NumAromaticRings', 'CeqC_DoubleBonds', 'Aldehyde_flag', 'Ketone_flag', 'CarboxAcid_flag',
         'Phenol_flag', 'Ether_flag', 'EsterLacton_flag', 'Alcohol_flag']
def groups(A):
    g = collections.defaultdict(list)
    for i in range(len(A)): g[np.ascontiguousarray(A[i]).tobytes()].append(i)
    return [v for v in g.values() if len(v) > 1], len(g)
gr, nd = groups(RB[RED12].to_numpy(float)); yc = list(RB.columns[5:10])
res['routeB_identical_rows'] = dict(distinct_vectors=nd, compounds=len(RB), compounds_sharing_a_vector=sum(len(v) for v in gr),
    groups=[dict(compounds=[RB['Monoterpenoid'][i] for i in v],
                 strains_with_different_classes=[c for c in yc if len({str(RB[c][i]).strip() for i in v}) > 1]) for v in gr])

# ---- benchmark, main analysis
bm = {}
for f in sorted(glob.glob(os.path.join(R, 'benchmark', 'datasets_k', 'DS*_X.csv'))):
    ds = os.path.basename(f)[:4]
    if int(ds[2:]) > 16: continue
    A = pd.read_csv(f, float_precision='round_trip').select_dtypes('number').to_numpy(float); y = pd.read_csv(f.replace('_X.csv', '_y.csv'))
    rc = y['relation'].astype(str).str.strip().isin(['>', '>=']); gr, _ = groups(A)
    bm[ds] = dict(n=int(len(y)), right_censored=int(rc.sum()), right_censored_percent=round(100 * rc.mean(), 1), right_censored_active=int(y['y'][rc].sum()),
                  identical_row_groups=[len(v) for v in gr], identical_rows_same_class=all(len(set(y['y'][v])) == 1 for v in gr))
res['benchmark'] = bm

# ---- peptides
rows = json.load(open(os.path.join(R, 'revision', 'lone', 'rows.json'))); P = pd.read_csv(os.path.join(R, 'revision', 'lone', 'desc_pool.csv'))
ids = [r['No'] for r in rows]; assert list(P['No'].astype(str)) == ids; X = P.iloc[:, 1:].to_numpy(float); cn = list(P.columns[1:])
def cls(T): return {k: np.array([1 if (float(r['MIC_' + k]) <= T[k] and r['censored_' + k] != 'yes') else 0 for r in rows]) for k in T}
YA, YC = cls({'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}), cls({'EC': 16, 'SA': 16, 'PA': 16, 'AB': 16})
mw = np.array([r['MW'] for r in rows]); tr = np.array([r['nres'] for r in rows]) < 8; sdp = X.std(axis=0)
i3 = [ids.index(k) for k in ('30', '31', '32')]; o3 = [i for i in range(len(ids)) if i not in i3]
res['peptides'] = dict(active_rule_A={k: int(v.sum()) for k, v in YA.items()}, active_rule_C={k: int(v.sum()) for k, v in YC.items()},
    P_aeruginosa_same_classes_under_both_rules=bool((YA['PA'] == YC['PA']).all()),
    truncated=[r['No'] for r, t in zip(rows, tr) if t], truncated_active_rule_A={k: int(v[tr].sum()) for k, v in YA.items()},
    truncated_non_active_against_every_strain=[r['No'] for i, r in enumerate(rows) if tr[i] and all(YA[k][i] == 0 for k in YA)],
    auc_of_molecular_weight_alone_rule_A={k: round(float(roc_auc_score(v, mw)), 3) for k, v in YA.items()},
    threshold_8_ug_per_mL_in_uM=[round(8000 / mw.max(), 1), round(8000 / mw.min(), 1)],
    columns_sd_below_1e_6={cn[j]: float(sdp[j]) for j in np.where(sdp < 1e-6)[0]},
    peptides_30_to_32_equal_to_median_of_the_other_34=[c for j, c in enumerate(cn) if len(set(X[i3, j])) == 1 and X[i3[0], j] == np.median(X[o3, j])
                                                    and c in ('Mordred::MAXsCH3', 'Mordred::MINsCH3', 'Mordred::MDEC-12', 'Mordred::MDEC-13')])
json.dump(res, open(OUT, 'w'), indent=1, ensure_ascii=False)
m = res['missing_monoterpenoids']; print('missing', m['missing_cells'], 'of', m['cells'], m['percent'], '%', m['imputed_columns_by_library'], len(m['compounds']), 'compounds')
print('replicates', {s: (v['active_mean_rule'], len(v['active_median_rule']), len(v['active_on_a_minority_of_replicates'])) for s, v in rep.items()}, res['replicates']['effective_endpoints_four_evaluable'])
print('Route B', nd, 'distinct of', len(RB), '|', res['routeB_identical_rows']['compounds_sharing_a_vector'], 'share a vector')
print('benchmark censored', {k: v['right_censored'] for k, v in bm.items() if v['right_censored']}, 'identical rows', {k: v['identical_row_groups'] for k, v in bm.items() if v['identical_row_groups']})
print('peptides', res['peptides']['auc_of_molecular_weight_alone_rule_A'], res['peptides']['truncated_non_active_against_every_strain'], res['peptides']['columns_sd_below_1e_6'])
