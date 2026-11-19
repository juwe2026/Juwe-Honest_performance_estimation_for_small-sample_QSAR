# -*- coding: utf-8 -*-
"""check_af2.py - compares every sheet of Additional file 2 with its source.

Result sheets are compared with the result objects of this archive (Additional file 10), the table sheets cell by
cell with Additional file 1 and, when given, with the main text. A printed number counts as matching its source when
it equals the source value rounded to the number of decimals printed. Nothing is recomputed from scratch except where
stated (null-distribution summaries, Clopper-Pearson intervals, effective number of endpoints, class labels from the
replicate MIC values), and those recomputations use only files of the archive.

Usage:
  python check_af2.py AF2.xlsx ARCHIVE_DIR AF1.docx TABLE_S10.json TIEBREAK.json [--base base.pkl] [--main MAIN.docx]

  ARCHIVE_DIR     unpacked code archive (code/, revision/, results/)
  TABLE_S10.json  output of revision/cp_borderline_all_rev2.py
  TIEBREAK.json   output of revision/jobB_tiebreak_0167.py
  base.pkl        rebuilt by revision/make_base.py, needed only for the descriptor names of Correlation_groups_182
  MAIN.docx       main text, needed only for Main_Table_1 to Main_Table_3

Exit code 0 when every check passes. Sheets whose source is not given are reported as not checked.

Numbering: Additional file 2 uses the numbering of the revised manuscript in citation order. The result objects of
the archive were written with the working numbers of the revision. Text taken from them is converted with
nummern_freigegeben.rn_old_text before it is compared (former main Tables 2 to 4 are Tables S16 to S18).
"""
import sys, os, re, json, glob, argparse, pickle
from decimal import Decimal
import numpy as np
import openpyxl, docx
from docx.oxml.ns import qn
from scipy.stats import beta
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nummern_freigegeben as NUMF

ap = argparse.ArgumentParser()
ap.add_argument('af2'); ap.add_argument('arch'); ap.add_argument('af1'); ap.add_argument('s10'); ap.add_argument('tiebreak')
ap.add_argument('--base'); ap.add_argument('--main')
A = ap.parse_args()
wb = openpyxl.load_workbook(A.af2)
RES = os.path.join(A.arch, 'results'); RV = os.path.join(RES, 'revision'); RP = os.path.join(RES, 'raw_permutations')
KEY = {'H. influenzae': 'HI', 'S. aureus': 'SA', 'S. pneumoniae': 'SPneu', 'S. pyogenes': 'SPyo', 'P. aeruginosa': 'PA'}
PKEY = {'E. coli ATCC 25922': 'EC', 'S. aureus ATCC 29213': 'SA', 'P. aeruginosa ATCC 27853': 'PA', 'A. baumannii ATCC 19606': 'AB'}
J = lambda *p: json.load(open(os.path.join(*p)))
REPORT, FAIL = [], []

def decimals(v):
    e = Decimal(repr(float(v))).normalize().as_tuple().exponent
    return max(0, -e)

def same(sheet_v, src_v):
    """sheet value equals the source rounded to the printed decimals (integers exactly)"""
    if sheet_v is None or src_v is None: return sheet_v is None and src_v is None
    if isinstance(sheet_v, str) or isinstance(src_v, str): return str(sheet_v) == str(src_v)
    if isinstance(sheet_v, bool) or isinstance(src_v, bool): return sheet_v == src_v
    if float(sheet_v) == float(src_v): return True
    nd = decimals(sheet_v)
    return abs(float(sheet_v) - float(src_v)) <= 0.5 * 10 ** (-nd) + 1e-12

def rev2_names(s):
    """Source text written before two approved renamings of 07.11.2026: the second data set is called 'Peptides'
    (item F01), and P. aeruginosa values are marked '(descriptive)' (item D2). Applied to source text only."""
    return s.replace('Octapeptides', 'Peptides').replace('(artefact)', '(descriptive)')

class Sheet:
    def __init__(self, name):
        self.name, self.n, self.bad = name, 0, []
    def eq(self, where, sheet_v, src_v):
        self.n += 1
        if isinstance(src_v, str): src_v = rev2_names(NUMF.rn_old_text(src_v))   # source text in working numbers
        if not same(sheet_v, src_v): self.bad.append('%s: sheet %r, source %r' % (where, sheet_v, src_v))
    def ok(self, where, cond, msg=''):
        self.n += 1
        if not cond: self.bad.append('%s %s' % (where, msg))
    def done(self, source):
        REPORT.append((self.name, self.n, len(self.bad), source))
        if self.bad: FAIL.append((self.name, self.bad[:12]))

def rows(name, start=2):
    ws = wb[name]; hdr = [c.value for c in ws[1]]
    return hdr, [dict(zip(hdr, r)) for r in ws.iter_rows(min_row=start, values_only=True)]

# ------------------------------------------------------------------ Seeds
S = Sheet('Seeds'); ws = wb['Seeds']
PYS = {}
for p in glob.glob(os.path.join(A.arch, 'code', '**', '*.py'), recursive=True) + glob.glob(os.path.join(A.arch, 'revision', '**', '*.py'), recursive=True):
    rel = os.path.relpath(p, A.arch); PYS.setdefault(os.path.basename(p), []).append(p); PYS[rel.replace('revision/', '', 1)] = [p]
for r in range(2, ws.max_row + 1):
    a, b, c, dv, e, f = [ws.cell(r, k).value for k in range(1, 7)]
    src = ''
    for tok in [t.strip() for t in re.split(r' / |, ', b)]:
        m = re.match(r'^([\w/]+\.py)', tok)
        if m:
            path = PYS.get(m.group(1)) or PYS.get(os.path.basename(m.group(1)))
            S.ok('row %d script %s' % (r, m.group(1)), bool(path), 'missing from the archive')
            if path: src += open(path[0]).read()
        m = re.match(r'^(\w+)\.(\w+)', tok)
        if m and not tok.endswith('.py') and not tok.startswith(('B_hits', 'E_misc', 'C_routeB')):
            mod = PYS.get(m.group(1) + '.py')
            S.ok('row %d function %s' % (r, tok.split(',')[0]), bool(mod) and ('def %s(' % m.group(2)) in open(mod[0]).read(), 'not found')
    if '(inherited)' not in str(dv) and 'as above' not in str(dv):
        for sd in re.findall(r'\d+', str(dv)):
            if re.search(r'\d+-\d+', str(dv)) and sd not in re.findall(r'\d+', str(dv))[:1]: continue
            S.ok('row %d seed %s' % (r, sd), re.search(r'(?<![\d.])%s(?![\d.])' % sd, src) is not None, 'not in the source code')
    elif 'as above' in str(dv):
        S.ok('row %d seed +100' % r, '+100' in src.replace(' ', ''), 'not in the source code')
    # files of column F
    base = None
    for tok in re.findall(r'[\w./*-]+\.(?:npz|json|npy|csv)|[\w./-]+/(?=[.\s,]|$)', f):
        if tok.startswith('results/'): base = os.path.dirname(tok); path = tok
        elif '/' in tok and not tok.startswith('code/'): path = os.path.join(base or 'results/revision', tok)
        elif tok.startswith('code/'): path = tok
        else: path = os.path.join(base, tok)
        hits = glob.glob(os.path.join(A.arch, path))
        S.ok('row %d file %s' % (r, path), len(hits) > 0, 'not in the archive')
        if isinstance(e, int) and path.endswith('.npz') and ('raw_permutations' in path or 'revision_nulls' in path or 'permB_25_cineole0' in path):
            allowed = {e} | ({500} if '500 permutations for the random forest' in f else set())
            lens = set()
            for h in hits:
                z = np.load(h); k = [x for x in z.files if x in ('q2s', 'aucs', 'null')] or [x for x in z.files if z[x].ndim == 1 and z[x].size > 10]
                if k: lens.add(int(z[k[0]].size))
            S.ok('row %d permutations %s' % (r, path), lens <= allowed and e in lens, 'length %s instead of %s' % (sorted(lens), sorted(allowed)))
S.done('source code and files of the archive')

# ------------------------------------------------------------------ Fallback_branch_audit
S = Sheet('Fallback_branch_audit'); FA = J(RV, 'fallback_audit.json')
hdr, rs = rows('Fallback_branch_audit'); S.ok('rows', len(rs) == len(FA), '%d versus %d' % (len(rs), len(FA)))
for x in rs:
    src = FA['%s | %s' % (x['Arm'], KEY[x['Endpoint']])]
    w = '%s %s' % (x['Arm'], x['Endpoint'])
    S.eq(w + ' Seed', x['Seed'], src['seed']); S.eq(w + ' n', x['Permutations_counted'], src['nperm'])
    S.eq(w + ' observed', x['Folds_taking_branch_observed'], '%d of 31' % src['folds_observed'])
    S.eq(w + ' permutiert %', x['Folds_taking_branch_permuted_pct'], src['folds_permuted_pct'])
    S.eq(w + ' Permutationen %', x['Permutations_with_any_such_fold_pct'], src['permutations_affected_pct'])
S.done('results/revision/fallback_audit.json')

# ------------------------------------------------------------------ Number_check
S = Sheet('Number_check'); NC = J(RV, 'numbercheck_rev.json')
hdr, rs = rows('Number_check'); S.ok('rows', len(rs) == NC['n'] == len(NC['rows']))
FULL_M = {v: k for k, v in KEY.items()}; FULL_P = {v: k for k, v in PKEY.items()}
asc = lambda t: t.replace('\u2212', '-') if isinstance(t, str) else t          # the sheets print minus signs as ASCII hyphens
for x, y in zip(rs, NC['rows']):
    w = '%s %s %s' % (x['Table'], x['Row'], x['Endpoint'])
    for a, b in (('Table', 'Table'), ('Row', 'Row'), ('Printed', 'Printed'), ('Source_value', 'Source_value'), ('Source_p', 'Source_p')): S.eq(w + ' ' + a, asc(x[a]), asc(y[b]))
    full = (FULL_P if y['Table'] == 'Table 5' else FULL_M).get(y['Endpoint'], y['Endpoint'])
    S.ok(w, x['Endpoint'] in (y['Endpoint'], full, full.split(' ATCC')[0]), 'endpoint %r versus %r' % (x['Endpoint'], y['Endpoint']))
    for a, b in (('Value_matches', 'Value_ok'), ('p_matches', 'p_ok'), ('Marker_matches', 'Marker_ok')): S.eq(w + ' ' + a, x[a], 'yes' if y[b] else 'no')
S.done('results/revision/numbercheck_rev.json (Tables 1 and 2, Tables S16 to S18, Quelle in Arbeitsnummern, umgerechnet)')

# ------------------------------------------------------------------ Null_distribution_audit (recomputed)
S = Sheet('Null_distribution_audit'); hdr, rs = rows('Null_distribution_audit')
for x in rs:
    f = os.path.join(RP, x['File']); S.ok(x['File'], os.path.exists(f), 'missing')
    if not os.path.exists(f): continue
    z = np.load(f); q, au = z['q2s'], z['aucs']; bl = np.abs(q) > 10
    for col, v in (('Permutations', len(q)), ('Permutations_with_abs_Q2_above_10', int(bl.sum())), ('Share_with_abs_Q2_above_10_pct', 100 * bl.mean()),
                   ('Q2_median', float(np.median(q))), ('Q2_mean_without_abs_Q2_above_10', float(q[~bl].mean())), ('AUC_mean', float(au.mean()))):
        S.eq('%s %s' % (x['File'], col), x[col], v)
S.done('results/raw_permutations/*.npz, nachgerechnet')

# ------------------------------------------------------------------ Monte Carlo sheets (Clopper-Pearson recomputed)
T10 = json.load(open(A.s10))
def cp(k, n):
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1); hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return lo, hi
def thr_in(lo, hi): return ' and '.join(t for t, v in (('0.05', 0.05), ('0.0167', 0.05 / 3)) if lo <= v <= hi)
S = Sheet('MonteCarlo_uncertainty'); hdr, rs = rows('MonteCarlo_uncertainty')
ill = [x for x in rs if x['Kind'] == 'illustration']; dat = [x for x in rs if x['Kind'] != 'illustration']
for x in ill:
    lo, hi = cp(x['k'], x['Permutations']); w = 'illustration n=%s k=%s' % (x['Permutations'], x['k'])
    S.eq(w + ' p', x['p_reported'], (x['k'] + 1) / (x['Permutations'] + 1)); S.eq(w + ' lo', x['CI95_low'], lo); S.eq(w + ' hi', x['CI95_high'], hi)
    S.eq(w + ' min p', x['Smallest_attainable_p'], 1 / (x['Permutations'] + 1))
S.ok('data rows', len(dat) == len(T10['table_S10']), '%d versus %d' % (len(dat), len(T10['table_S10'])))
for x in dat:
    lo, hi = cp(x['k'], x['Permutations']); w = x['Kind']
    S.eq(w + ' p', x['p_reported'], (x['k'] + 1) / (x['Permutations'] + 1)); S.eq(w + ' lo', x['CI95_low'], lo); S.eq(w + ' hi', x['CI95_high'], hi)
S.done('%s, Intervalle nachgerechnet' % os.path.basename(A.s10))
S = Sheet('MonteCarlo_all_checked'); hdr, rs = rows('MonteCarlo_all_checked')
S.ok('rows', len(rs) == T10['n_checked'] == len(T10['items']))
src_k = sorted((rev2_names(i['Data']), i['Model'].replace('Route A1u hit count, tie-break largest variance, as published', 'Route A1u hit count, tie-break largest variance'), i['Endpoint'], i['Metric'], i['k'], i['n']) for i in T10['items'])
sh_k = sorted((x['Data_set'], x['Model'], x['Endpoint'], x['Metric'], x['k'], x['Permutations']) for x in rs)
S.ok('Inhalt', src_k == sh_k, 'k and n differ from the source')
for x in rs:
    lo, hi = cp(x['k'], x['Permutations']); w = '%s %s %s %s' % (x['Source'], x['Model'], x['Endpoint'], x['Metric'])
    S.eq(w + ' p', x['p_reported'], (x['k'] + 1) / (x['Permutations'] + 1)); S.eq(w + ' lo', x['CI95_low'], lo); S.eq(w + ' hi', x['CI95_high'], hi)
    S.eq(w + ' threshold', x['Interval_contains'] or '', thr_in(lo, hi))
S.done('%s, intervals and thresholds recomputed' % os.path.basename(A.s10))

# ------------------------------------------------------------------ Table1_all_arms
S = Sheet('Table1_all_arms'); hdr, rs = rows('Table1_all_arms')
IF, HI2, NF, RF, EM, M2, B25 = (J(RV, 'I_fallback_seeds.json'), J(RV, 'A2_HI_seed42.json'), J(RV, 'A_nofilter.json'), J(RV, 'D_rf.json'),
                                J(RV, 'E_misc.json'), J(RV, 'morgan2_fixed.json'), J(RV, 'B25_cineole0_seed701.json'))
def src_t1(model, k):
    if model.startswith('Route A1 (') or model.startswith('Route A1u (') or model.startswith('Route A2 ('):
        rt = model.split()[1]; s = HI2 if (rt == 'A2' and k == 'HI') else IF[k][rt]
        return dict(q2=s['q2'], auc=s['auc'], p_q2=s['p_q2'], p_auc=s['p_auc'], n=s['nperm'], ci=EM['nested_obs'][k][rt]['ci95'])
    m = {'No selection, all 2,208': NF['all2208'], 'No filter, 182 label-free representatives': NF['reps182_centrality'],
         'No filter, 182 Route A1 representatives (nested)': NF['reps182_effectsize_nested']}
    if model in m: s = m[model][k]; return dict(q2=s['q2'], auc=s['auc'], p_q2=s['p_q2'], p_auc=s['p_auc'], n=s['nperm'])
    m = {'Random forest, out-of-bag, all 2,208': 'RouteA_all2208', 'Random forest, out-of-bag, 182 representatives': 'RouteA_182_centrality',
         'Random forest, out-of-bag, Route B 12': 'RouteB_12'}
    if model in m: s = RF[m[model]][k]; return dict(q2=s['oob_q2'], auc=s['oob_auc'], p_q2=s['p_q2'], p_auc=s['p_auc'], n=s['nperm'])
    if model == 'Morgan2 fingerprints, no selection': s = M2[k]['no_selection']; return dict(q2=s['q2'], auc=s['auc'], p_q2=s['p_q2'], p_auc=s['p_auc'], n=2000)
    if model == 'Morgan2 fingerprints, nested filter': s = M2[k]['nested_filter']; return dict(q2=s['q2'], auc=s['auc'], p_q2=s['p_q2'], p_auc=s['p_auc'], n=2000)
    if model.startswith('Route B, 25 descriptors'): s = B25[k]; return dict(q2=s['q2'], auc=s['auc'], p_q2=s['p_q2'], p_auc=s['p_auc'], n=s['n_perm'])
    raise KeyError(model)
for x in rs:
    k = KEY[x['Strain']]; s = src_t1(x['Model'], k); w = '%s %s' % (x['Strain'], x['Model'])
    S.eq(w + ' Q2', x['Q2'], s['q2']); S.eq(w + ' AUC', x['AUC'], s['auc']); S.eq(w + ' n', x['Permutations'], s['n'])
    if k == 'PA' and x['p_Q2'] is None: S.ok(w, x['Descriptive_only'] == 'yes', 'p missing without descriptive-only flag')
    else: S.eq(w + ' p_Q2', x['p_Q2'], s['p_q2']); S.eq(w + ' p_AUC', x['p_AUC'], s['p_auc'])
    if 'ci' in s: S.eq(w + ' CI lo', x['Q2_CI95_low'], s['ci'][0]); S.eq(w + ' CI hi', x['Q2_CI95_high'], s['ci'][1])
    S.ok(w + ' Descriptive_only', (x['Descriptive_only'] == 'yes') == (k == 'PA'))
S.done('I_fallback_seeds, A2_HI_seed42, A_nofilter, D_rf, E_misc, morgan2_fixed, B25_cineole0_seed701')

# ------------------------------------------------------------------ HitCounts_all_data_sets
S = Sheet('HitCounts_all_data_sets'); hdr, rs = rows('HitCounts_all_data_sets')
BH, GH, LA, LC = J(RV, 'B_hits.json'), J(RV, 'G_gaucher_hits.json'), J(RV, 'lone', 'res_A.json'), J(RV, 'lone', 'res_C.json')
for x in rs:
    rt = x['Route'].replace('Route ', '')
    if x['Data_set'].startswith('Monoterpenoids'): s = dict(BH['calibration'][KEY[x['Endpoint']]][{'A1': 'A1s'}.get(rt, rt)]); s.setdefault('n_tests', {'A2': 2208}.get(rt, 182))
    elif x['Data_set'].startswith('Gaucher'): s = dict(GH[rt]); s.setdefault('n_tests', GH['n_var'] if rt == 'A2' else GH['n_groups'])
    else: s = (LA if 'main' in x['Data_set'] else LC)[PKEY[x['Endpoint']]]['hits'][rt]
    w = '%s %s %s' % (x['Data_set'][:12], x['Endpoint'], rt)
    for col, key in (('n_tests', 'n_tests'), ('Observed_hits', 'obs'), ('Null_mean', 'null_mean'), ('Null_SD', 'null_sd'), ('Null_hit_rate', 'rate'),
                     ('Inflation_vs_nominal_5pct', 'infl_nominal'), ('Inflation_vs_own_rate', 'infl_own'), ('p', 'p')):
        S.eq(w + ' ' + col, x[col], s[key])
if 'p_B100000' in hdr:                     # columns from revision2/analysis/hits_extended.py (B = 100,000)
    def HX(part):
        f = sorted(glob.glob(os.path.join(A.arch, 'revision2', 'results', 'hits_extended_%s*.json' % part)))
        S.ok('hits_extended_%s result file' % part, len(f) == 1, '%d files' % len(f)); return json.load(open(f[-1]))['results'] if f else {}
    HM, HP, HG = HX('monoterpenoids'), HX('peptides'), HX('gaucher')
    KMX = {'H. influenzae': 'H.I.', 'S. aureus': 'S.A.', 'S. pneumoniae': 'S.Pneu', 'S. pyogenes': 'S.Pyo', 'P. aeruginosa': 'P.A.'}
    for x in rs:
        rt = x['Route'].replace('Route ', ''); w = '%s %s %s' % (x['Data_set'][:12], x['Endpoint'], rt)
        if x['Data_set'].startswith('Monoterpenoids'): v = HM[KMX[x['Endpoint']]][rt]
        elif x['Data_set'].startswith('Gaucher'): v = HG[rt]
        else: v = HP['A' if 'main' in x['Data_set'] else 'C'][PKEY[x['Endpoint']]][rt]
        S.eq(w + ' p (B1000) from the extension', x['p'], v['B1000']['p'])
        for col, val in (('k_B1000', v['B1000']['k']), ('CI95_low_k_over_B_B1000', round(v['B1000']['ci95_k_over_B'][0], 5)),
                         ('CI95_high_k_over_B_B1000', round(v['B1000']['ci95_k_over_B'][1], 5)), ('p_B100000', round(v['B100000']['p'], 6)),
                         ('k_B100000', v['B100000']['k']), ('CI95_low_k_over_B_B100000', round(v['B100000']['ci95_k_over_B'][0], 6)),
                         ('CI95_high_k_over_B_B100000', round(v['B100000']['ci95_k_over_B'][1], 6))):
            S.eq(w + ' ' + col, x[col], val)
        k = v['B1000']['k']; lo = 0 if k == 0 else beta.ppf(0.025, k, 1000 - k + 1); hi = beta.ppf(0.975, k + 1, 1000 - k)
        S.ok(w + ' interval B1000 recomputed', abs(lo - x['CI95_low_k_over_B_B1000']) < 1e-5 and abs(hi - x['CI95_high_k_over_B_B1000']) < 1e-5)
        S.ok(w + ' p B1000 = (k + 1)/1001', abs((k + 1) / 1001 - x['p']) < 6e-5)
S.done('B_hits, G_gaucher_hits, lone/res_A, lone/res_C, revision2/results/hits_extended_*.json')

# ------------------------------------------------------------------ SecondDataSet_models
S = Sheet('SecondDataSet_models'); hdr, rs = rows('SecondDataSet_models')
LAF, LCF = J(RV, 'lone', 'res_A_fix.json'), J(RV, 'lone', 'res_C_fix.json')
MOD = {'Route A1 (nested)': 'A1', 'Route A1u (nested)': 'A1u', 'Route A2 (nested)': 'A2', 'No selection, all 1,275': 'noSelection_all',
       'No filter, 50 label-free representatives': 'noFilter_reps', 'Random forest, out-of-bag': 'randomForest_oob'}
for x in rs:
    main = x['Definition'].startswith('main'); res, fix = (LA, LAF) if main else (LC, LCF)
    k = PKEY[x['Strain']]; m = MOD[x['Model']]; s = (fix if m in ('A1', 'A1u', 'A2') else res)[k][m]
    w = '%s %s %s' % ('A' if main else 'C', k, m)
    S.eq(w + ' Schwelle', x['Threshold_ug_mL'], res['thresholds'][k]); S.eq(w + ' n_active', x['n_active'], res['n_active'][k]); S.eq(w + ' n', x['n_total'], 37)
    for col, key in (('Q2', 'q2'), ('p_Q2', 'p_q2'), ('AUC', 'auc'), ('p_AUC', 'p_auc'), ('MCC', 'mcc'), ('TP', 'TP'), ('FN', 'FN'), ('FP', 'FP'), ('TN', 'TN')):
        S.eq(w + ' ' + col, x[col], s[key])
    S.eq(w + ' n_perm', x['Permutations'], s.get('nperm', res['nperm']))
S.done('lone/res_A(_fix), lone/res_C(_fix)')

# ------------------------------------------------------------------ Correlation_groups_182
S = Sheet('Correlation_groups_182'); hdr, rs = rows('Correlation_groups_182')
G = pickle.load(open(os.path.join(RP, 'alt_groups.pkl'), 'rb')); U = set(int(i) for i in pickle.load(open(os.path.join(RP, 'unsup_reps.pkl'), 'rb')))
S.ok('groups', len(G) == 182 and sum(len(g) for g in G) == 2208 == len(rs))
if A.base:
    DESC = pickle.load(open(A.base, 'rb'))['D']['HI']['DESC']
    byg = {}
    for x in rs: byg.setdefault(x['group'], []).append(x)
    S.ok('number of groups in the sheet', sorted(byg) == list(range(1, 183)))
    for gi, g in enumerate(G, 1):
        xs = byg.get(gi, []); names = [DESC[i] for i in g]; rep = [DESC[i] for i in g if i in U]
        S.ok('group %d members' % gi, sorted(x['member'] for x in xs) == sorted(names))
        S.ok('group %d size' % gi, all(x['size'] == len(g) for x in xs)); S.ok('group %d hub rank' % gi, all(x['hub_order'] == gi for x in xs))
        S.ok('group %d representative' % gi, len(rep) == 1 and all(x['representative_A1u_centrality'] == rep[0] for x in xs), str(rep[:2]))
    S.done('raw_permutations/alt_groups.pkl, unsup_reps.pkl, names from base.pkl')
else:
    S.done('raw_permutations/alt_groups.pkl (without base.pkl only number and size)')

# ------------------------------------------------------------------ TieBreak sheets
TB = json.load(open(A.tiebreak))
S = Sheet('TieBreak_sensitivity'); hdr, rs = rows('TieBreak_sensitivity')
LBL = {'largest variance (used in this work)': 'variance (as published)', 'first index': 'first index'}
for x in rs:
    k = KEY[x['Endpoint']]; w = '%s %s' % (x['Tie_break_rule'], k)
    if x['Tie_break_rule'] in LBL:
        s = BH['tiebreak_fixed'][LBL[x['Tie_break_rule']]]
        S.eq(w + ' groups', x['n_groups'], s['n_groups'])
        for col, key in (('Observed_hits', 'obs'), ('Null_mean', 'null_mean'), ('Null_SD', 'null_sd'), ('p', 'p')): S.eq(w + ' ' + col, x[col], s[k][key])
    else:
        s = TB['summary'][k]; ng = [TB['per_rule'][r]['n_groups'] for r in TB['per_rule']]
        S.eq(w + ' groups', x['n_groups'], '%d-%d (median %d)' % (min(ng), max(ng), np.median(ng)))
        obs = [TB['per_rule'][r][k]['obs'] for r in TB['per_rule']]
        S.eq(w + ' hits', x['Observed_hits'], '%d-%d (median %d)' % (min(obs), max(obs), np.median(obs)))
        S.eq(w + ' p', x['p'], '%.4f-%.4f (median %.4f), p<=0.05 in %.0f%% of rules, p<=0.0167 in %.0f%%' % (s['p_min'], s['p_max'], s['p_median'], 100 * s['frac_p_le_05'], 100 * s['frac_p_le_0167']))
S.done('B_hits.tiebreak_fixed, %s' % os.path.basename(A.tiebreak))
S = Sheet('TieBreak_random_rules'); hdr, rs = rows('TieBreak_random_rules')
S.ok('rows', len(rs) == 5 * len(TB['per_rule']))
for x in rs:
    s = TB['per_rule'][x['Rule']]; k = KEY[x['Endpoint']]; w = '%s %s' % (x['Rule'], k)
    S.eq(w + ' groups', x['n_groups'], s['n_groups'])
    for col, key in (('Observed_hits', 'obs'), ('Null_mean', 'null_mean'), ('Null_SD', 'null_sd'), ('p', 'p')): S.eq(w + ' ' + col, x[col], s[k][key])
for k in KEY.values():                                   # the summary is consistent with the 100 rules
    ps = np.array([TB['per_rule'][r][k]['p'] for r in TB['per_rule']])
    S.eq(k + ' Median', round(float(np.median(ps)), 4), round(TB['summary'][k]['p_median'], 4))
    S.eq(k + ' share 0.0167', float(np.mean(ps <= 0.05 / 3 + 1e-12)), TB['summary'][k]['frac_p_le_0167'])
S.done(os.path.basename(A.tiebreak))

# ------------------------------------------------------------------ Activity_vectors_31, Endpoint_correlation, MIC_replicates_31
EC = J(RV, 'endpoint_corr.json'); mono = EC['mono']
S = Sheet('Activity_vectors_31'); hdr, rs = rows('Activity_vectors_31')
S.ok('compounds', [x['Compound'] for x in rs] == mono['compounds'])
for i, x in enumerate(rs):
    for nm, k in KEY.items(): S.eq('%s %s' % (x['Compound'], k), x[nm], mono['Y'][k][i])
S.done('results/revision/endpoint_corr.json')
S = Sheet('Endpoint_correlation'); ws = wb['Endpoint_correlation']; hdr = [c.value for c in ws[1]]
COLS = {h.replace('phi vs ', ''): i + 1 for i, h in enumerate(hdr) if h and h.startswith('phi vs ')}
SHORT = {'HI': 'H. influenzae', 'SA': 'S. aureus', 'SPneu': 'S. pneumoniae', 'SPyo': 'S. pyogenes', 'PA': 'P. aeruginosa', 'EC': 'E. coli', 'AB': 'A. baumannii'}
L4 = np.array([mono['Y'][k] for k in ('HI', 'SA', 'SPneu', 'SPyo')], float); C4 = np.corrcoef(L4); e4 = np.sort(np.linalg.eigvalsh(C4))[::-1]
four = dict(strains=['HI', 'SA', 'SPneu', 'SPyo'], phi=C4.tolist(), eigen=e4.tolist(), li_ji=float(sum((e >= 1) + (e - np.floor(e)) for e in e4)),
            nyholt=float(1 + 3 * (1 - np.var(e4, ddof=1) / 4)))
BLK = [('31 monoterpenoids (mean', mono), ('31 monoterpenoids, four', four), ('37 peptides, rule A', EC['lone_rule_A']), ('37 peptides, rule C', EC['lone_rule_C'])]
for prefix, src in BLK:
    rr = [r for r in range(2, ws.max_row + 1) if ws.cell(r, 1).value.startswith(prefix)]
    st = src['strains']; S.ok(prefix + ' rows', len(rr) == len(st) + 2)
    for j, k in enumerate(st):
        r = rr[j]; S.ok('%s %s' % (prefix, k), ws.cell(r, 2).value == SHORT[k])
        for m, k2 in enumerate(st): S.eq('%s phi %s-%s' % (prefix, k, k2), ws.cell(r, COLS[SHORT[k2]]).value, src['phi'][j][m])
    r = rr[len(st)]; vals = [ws.cell(r, c).value for c in range(3, ws.max_column + 1) if ws.cell(r, c).value is not None]
    S.ok(prefix + ' eigenvalues descending', vals == sorted(vals, reverse=True))
    for v, s in zip(vals, src['eigen']): S.eq(prefix + ' Eigenwert', v, s)
    S.eq(prefix + ' Meff', ws.cell(rr[-1], 2).value, 'effective number of endpoints: Li and Ji %.2f, Nyholt/Cheverud %.2f' % (src['li_ji'], src['nyholt']))
S.done('endpoint_corr.json, block of the four endpoints recomputed')
S = Sheet('MIC_replicates_31'); hdr, rs = rows('MIC_replicates_31')
liq = {}
for x in rs:
    if x['Phase'].startswith('liquid'): liq.setdefault((x['Compound'], x['Strain']), []).append(x['MIC_ug_mL'])
S.ok('pairs liquid phase', len(liq) == 155, str(len(liq)))
for (cpd, st), v in liq.items():
    i = mono['compounds'].index(cpd); S.eq('%s %s class' % (cpd, KEY[st]), int(np.mean(v) < 2048), mono['Y'][KEY[st]][i])
nrep = [len(v) for v in liq.values()]
S.ok('6 to 12 values per pair', min(nrep) == 6 and max(nrep) == 12, '%d-%d' % (min(nrep), max(nrep)))
S.ok('9 values for 135 of 155 pairs', sum(n == 9 for n in nrep) == 135, str(sum(n == 9 for n in nrep)))
S.done('classes recomputed from the single values against endpoint_corr.json (the laboratory workbook is not part of the archive)')

# ------------------------------------------------------------------ table sheets against Additional file 1 and the main text
norm = lambda s: re.sub(r'\s+', ' ', str(s if s is not None else '')).replace('−', '-').replace('–', '-').strip()
def doc_tables(path, pat):
    d = docx.Document(path); heads = {}; last = None; ti = 0
    for el in d.element.body.iterchildren():
        if el.tag == qn('w:p'):
            m = re.match(pat, docx.text.paragraph.Paragraph(el, d).text)
            if m: last = m.group(1)
        elif el.tag == qn('w:tbl'): heads.setdefault(last, []).append(d.tables[ti]); ti += 1
    return heads
def cmp_table(S, ws, tab, transposed=False):
    # transposed=True: the rows of the Word table are the columns of the sheet
    for i, row in enumerate(tab.rows, 1):
        seen = []
        for c in row.cells:
            if c._tc in seen: continue
            seen.append(c._tc); j = len(seen)
            r, k = (j, i) if transposed else (i, j)
            v = ws.cell(r, k).value
            if r == 1 and k == 1 and v == 'col0': continue
            S.ok('%s row %d col %d' % (ws.title, r, k), norm(v) == norm(c.text), '%r versus %r' % (str(v)[:40], c.text[:40]))
    if transposed:
        S.ok(ws.title + ' number of rows', ws.max_row == len(tab.columns) and ws.max_column == len(tab.rows),
             '%d x %d versus %d x %d' % (ws.max_row, ws.max_column, len(tab.columns), len(tab.rows)))
    else:
        S.ok(ws.title + ' number of rows', ws.max_row == len(tab.rows), '%d versus %d' % (ws.max_row, len(tab.rows)))
H1 = doc_tables(A.af1, r'(Table S\d+)[: ]')
S = Sheet('Table_S1 bis Table_S26')
names = sorted([n for n in wb.sheetnames if re.match(r'Table_S\d+(_part2)?$', n)], key=lambda n: (int(re.findall(r'\d+', n)[0]), n))
for n in names:
    k = 'Table S' + re.findall(r'\d+', n)[0]
    if k == 'Table S13': continue
    tabs = H1.get(k, []); j = 1 if n.endswith('_part2') else 0
    S.ok(n, len(tabs) > j, 'table missing in Additional file 1')
    if len(tabs) > j: cmp_table(S, wb[n], tabs[j])
S.ok('every table of Additional file 1 has a sheet', all(('Table_S%s' % k[7:]) in wb.sheetnames for k in H1 if k))
S.done('Additional file 1 (%s), cell by cell' % os.path.basename(A.af1))
S = Sheet('Table_S13'); hdr, rs = rows('Table_S13'); PB = J(RV, 'paired_bootstrap_Q2.json')
S.ok('rows', len(rs) == PB['n_comparisons'] == len(PB['rows']), '%d versus %d' % (len(rs), PB['n_comparisons']))
DS = {'ds1': 'Monoterpenoids', 'ds2': 'Peptides'}
yn = lambda b: 'yes' if b else 'no'
for x, y in zip(rs, PB['rows']):
    w = '%s %s vs %s %s' % (x['Data_set'][:5], x['Model_1'], x['Model_2'], x['Endpoint'])
    S.ok(w + ' data set', x['Data_set'].startswith(DS.get(y['data_set'], y['data_set'])[:5]), '%r versus %r' % (x['Data_set'], y['data_set']))
    for a, b in (('Model_1', 'model_1'), ('Model_2', 'model_2'), ('Endpoint', 'endpoint'), ('Q2_model_1', 'q2_model_1'), ('Q2_model_2', 'q2_model_2'),
                 ('Difference', 'difference'), ('CI95_lower', 'ci95_lower'), ('CI95_upper', 'ci95_upper'), ('p_two_sided_bootstrap', 'p_two_sided')):
        S.eq(w + ' ' + a, x[a], y[b])
    for a, b in (('CI_excludes_zero', 'ci_excludes_zero'), ('Holm_5pct', 'holm_5pct'), ('Benjamini_Hochberg_5pct', 'benjamini_hochberg_5pct')): S.eq(w + ' ' + a, x[a], yn(y[b]))
S.ok('Holm / BH / CI-Zahlen', (sum(x['Holm_5pct'] == 'yes' for x in rs), sum(x['Benjamini_Hochberg_5pct'] == 'yes' for x in rs), sum(x['CI_excludes_zero'] == 'yes' for x in rs)) == (PB['n_holm'], PB['n_bh'], PB['n_ci_excludes_zero']))
S.done('results/revision/paired_bootstrap_Q2.json, row by row')
if A.main:
    H0 = doc_tables(A.main, r'(Table \d+)[: ]')
    S = Sheet('Main_Table_*')
    # Table 3 of the main text shows the hypotheses as columns. Its sheet keeps one row per hypothesis.
    TRANSPOSED_IN_MAIN_TEXT = {'Main_Table_3'}
    for n in [n for n in wb.sheetnames if n.startswith('Main_Table_')]:
        k = 'Table ' + n.split('_')[-1]; S.ok(n, k in H0, 'table missing in the main text')
        if k in H0: cmp_table(S, wb[n], H0[k][0], transposed=n in TRANSPOSED_IN_MAIN_TEXT)
    S.done('main text (%s), cell by cell' % os.path.basename(A.main))
else:
    REPORT.append(('Main_Table_*', 0, 0, 'NOT CHECKED: main text not given (at assembly)'))

# ------------------------------------------------------------------ whole workbook
S = Sheet('ganze Mappe')
known = {'Info', 'Seeds', 'Fallback_branch_audit', 'Number_check', 'Null_distribution_audit', 'MonteCarlo_uncertainty', 'MonteCarlo_all_checked',
         'Table1_all_arms', 'HitCounts_all_data_sets', 'SecondDataSet_models', 'Correlation_groups_182', 'TieBreak_sensitivity',
         'TieBreak_random_rules', 'Activity_vectors_31', 'Endpoint_correlation', 'MIC_replicates_31'}
for n in wb.sheetnames:
    S.ok('sheet %s checked' % n, n in known or n.startswith(('Table_S', 'Main_Table_')), 'no check step')
HIST = re.compile(r'\b(submitted|reviewer|corrected|reconstructed|as published|recomputation|as deposited)\b', re.I)
GERMAN = re.compile(r'\b(Datei|anteil|mittel|ohne|und|nicht|Zeile)\b')
for n in wb.sheetnames:
    for row in wb[n].iter_rows():
        for c in row:
            if isinstance(c.value, str):
                if ';' in c.value and not (n == 'Table_S7' and c.value.startswith('[')): S.ok('%s!%s semicolon' % (n, c.coordinate), False, c.value[:60])
                if HIST.search(c.value) and n not in ('Table_S21',): S.ok('%s!%s Geschichte' % (n, c.coordinate), False, c.value[:60])
                if GERMAN.search(c.value): S.ok('%s!%s Deutsch' % (n, c.coordinate), False, c.value[:60])
info = {r[0]: r[1] for r in wb['Info'].iter_rows(min_row=2, values_only=True)}
S.ok('Info nennt check_af2.py', 'check_af2.py' in (info.get('Provenance') or ''))
S.done('sheet list, semicolons, references to earlier versions, German words')

# ------------------------------------------------------------------ report
w = max(len(r[0]) for r in REPORT)
for name, n, nb, src in REPORT:
    print('%-*s  %6d checks  %4d deviations   %s' % (w, name, n, nb, src))
for name, bad in FAIL:
    print('\n%s:' % name)
    for b in bad: print('   ' + b)
tot = sum(r[1] for r in REPORT); nb = sum(r[2] for r in REPORT)
print('\nTotal %d checks, %d deviations' % (tot, nb))
sys.exit(1 if nb else 0)
