# -*- coding: utf-8 -*-
"""Table S10, full scope (decision EW 06.11.2026, option A): Clopper-Pearson intervals (95 %) for the Monte Carlo
uncertainty of every permutation p-value reported in the main text and in Additional file 1 for the three data sets
of this study (monoterpenoids, octapeptides, Gaucher). Thresholds 0.05 and 0.05/3 = 0.0167.

Sources
  (1) Tables 1 to 5 (main text and the former Tables 2 to 4): the items of cp_borderline_rev2.py, exact p from JSON.
  (2) Additional file 1:
      Table S1  Route B, grouping recomputed in every fold, Q² 1 LV and 3 LV, 2,000 permutations
      Table S2  Routes A1 and A1u, groups fixed and rebuilt, nested Q², 2,000 (independent draw from Table 1)
      Table S3  median activity rule (H. influenzae, S. pyogenes), Q² and AUC, 2,000
      Table S5  tie-break rules of the grouping, Route A1u hit count, 1,000 per rule (independent nulls)
      Table S12 univariate hit counts, Routes A2, A1, A1u, 1,000
      Table S13 random forest (500) on the 182 representatives, Route B and Morgan bits, PLS-DA on Morgan bits (2,000)
      Table S20 Gaucher data, nested Q² and AUC, 1 and 2 LV, 2,000
      new section, Route A1 and A1u without the 3D descriptor blocks (F_no3d_fixed.json), 2,000
      p-values given only in the text of Additional file 1: Route B without the phenol indicator (11 descriptors,
      Q² for H. influenzae, S. aureus and S. pneumoniae, C_routeB.json, 2,000) and Route A2 restricted to its top-15
      descriptors (Q² for S. pneumoniae and S. pyogenes, cache_topk_seeds/*_A2_top15.json, 2,000)
  Values that repeat an entry of Tables 1 to 5 (Table S1 fixed rows, Table S3 mean-rule columns, Table S8 Route B
  3 LV, Table S13 all-descriptor forest) are checked for identity and counted once.
  P. aeruginosa of the monoterpenoids (three actives) is descriptive only and is excluded, as in Table 1.
  Apparent p-values reported as <= 0.001 (k = 0 of 1,000) cannot be borderline and are counted.
  Not included: the 100 random tie-break rules of Table S5 (reported as a range and as shares, see the note in the
  output), the hit counts of the octapeptides and the Gaucher data (not reported as p-values in the main text or
  Additional file 1), and the benchmark, whose hypothesis tests use no permutation p-values.

Usage: python revision2/analysis/cp_borderline_all_rev2.py ARCHIVE_ROOT AF2.xlsx OUT.json
"""
import json, os, re, sys, runpy, io, contextlib
import numpy as np
import openpyxl
from scipy.stats import beta

ROOT, AF2, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
O = os.path.join(ROOT, 'results', 'revision') + os.sep
THR = (0.05, 0.05 / 3)
J = lambda f: json.load(open(f))

def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)

def k_exact(p, n):
    k = p * (n + 1) - 1
    assert abs(k - round(k)) < 1e-6, (p, n)
    return int(round(k))

def k_rounded(s, n):
    """k from a p-value printed with d decimals: the unique k with (k+1)/(n+1) rounding to the printed value."""
    p = float(s); d = len(s.split('.')[1])
    cands = [k for k in range(n + 1) if abs(round((k + 1) / (n + 1), d) - p) < 10 ** -(d + 3)]
    assert len(cands) == 1, (s, n, cands)
    return cands[0]

# ---------------------------------------------------------------- (1) Tables 1 to 5
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(os.path.join(ROOT, 'revision', 'cp_borderline_rev2.py'),
                       init_globals=None, run_name='not_main') if False else None
sys_argv = sys.argv
sys.argv = ['x', ROOT, os.devnull]
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(os.path.join(ROOT, 'revision', 'cp_borderline_rev2.py'))
sys.argv = sys_argv
ITEMS = []
for it in g['ITEMS']:
    ds = 'Octapeptides' if it['Source'] == 'Table 5' else 'Monoterpenoids'
    ITEMS.append(dict(Data=ds, Source=it['Source'], Model=it['Model'], Endpoint=it['Endpoint'], Metric=it['Metric'],
                      p=it['p'], n=it['n'], k=k_exact(it['p'], it['n'])))
# apparent p-values of the former Tables 3 and 4 (<= 0.001, 1,000 permutations, k = 0)
for src, lab in (('Table 3', 'Route A1 top-15 by effect size, apparent'), ('Table 4', 'Route A1 top-15 by VIP, apparent')):
    for e in ('H. influenzae', 'S. aureus', 'S. pneumoniae', 'S. pyogenes'):
        for m in ('Q²', 'AUC'):
            ITEMS.append(dict(Data='Monoterpenoids', Source=src, Model=lab, Endpoint=e, Metric=m, p=1 / 1001, n=1000, k=0))

# ---------------------------------------------------------------- (2) Additional file 1
wb = openpyxl.load_workbook(AF2, data_only=True)
def rows(sheet):
    ws = wb[sheet]; it = ws.iter_rows(values_only=True); hdr = [str(c) if c is not None else '' for c in next(it)]
    return hdr, [[str(c) if c is not None else '' for c in r] for r in it]
PAT = re.compile(r'\((\d\.\d+)\)')
EP = {'H. influenzae': 'H. influenzae', 'S. aureus': 'S. aureus', 'S. pneumoniae': 'S. pneumoniae', 'S. pyogenes': 'S. pyogenes'}

def grid(sheet, src, ds, n, rowmap, colmap):
    """rowmap: row label prefix -> (model, metric). colmap: column header -> endpoint (None = skip)."""
    hdr, rs = rows(sheet); out = []
    for r in rs:
        hit = [v for pre, v in rowmap.items() if r[0].startswith(pre)]
        if not hit: continue
        model, metric = hit[0]
        for j, h in enumerate(hdr[1:], 1):
            e = colmap.get(h)
            if e is None or not r[j]: continue
            m = PAT.search(r[j])
            if not m: continue
            s = m.group(1)
            out.append(dict(Data=ds, Source=src, Model=model, Endpoint=e, Metric=metric, p=float(s), n=n, k=k_rounded(s, n)))
    return out

MONO = {h: h for h in EP}
# Table S1: Route B, grouping recomputed in every fold
ITEMS += grid('Table_S1', 'Table S1', 'Monoterpenoids', 2000,
              {'Q², 1 LV, grouping recomputed': ('Route B, 12 descriptors, 1 LV, grouping recomputed in every fold', 'Q²'),
               'Q², 3 LV, grouping recomputed': ('Route B, 12 descriptors, 3 LV, grouping recomputed in every fold', 'Q²')}, MONO)
# Table S2: Routes A1 and A1u, fixed and rebuilt groups (independent permutation draw)
ITEMS += grid('Table_S2', 'Table S2', 'Monoterpenoids', 2000,
              {'Route A1, nested Q², groups fixed': ('Route A1, groups fixed (paired draw)', 'Q²'),
               'Route A1, nested Q², groups rebuilt': ('Route A1, groups rebuilt in every fold', 'Q²'),
               'Route A1u, nested Q², groups fixed': ('Route A1u, groups fixed (paired draw)', 'Q²'),
               'Route A1u, nested Q², groups rebuilt': ('Route A1u, groups rebuilt in every fold', 'Q²')}, MONO)
# Table S3: median rule only (the mean-rule columns repeat Table 1 and are checked below)
S3rows = {'Route A1, nested': 'Route A1', 'Route A1u, nested': 'Route A1u', 'Route A2, nested': 'Route A2',
          'No selection, all 2,208 descriptors': 'No selection, all 2,208 descriptors',
          'Route B, 12 descriptors, 1 LV': 'Route B, 12 descriptors, 1 LV', 'Route B, 12 descriptors, 3 LV': 'Route B, 12 descriptors, 3 LV'}
hdr, rs = rows('Table_S3'); metric = None; s3_mean = []
for r in rs:
    if r[0].startswith('(a)'): metric = 'Q²'; continue
    if r[0].startswith('(b)'): metric = 'AUC'; continue
    if r[0] not in S3rows or metric is None: continue
    for j, h in enumerate(hdr[1:], 1):
        m = PAT.search(r[j]); e = h.split(',')[0]
        if not m: continue
        s = m.group(1)
        rec = dict(Data='Monoterpenoids', Endpoint=e, Metric=metric, p=float(s), n=2000, k=k_rounded(s, 2000))
        if 'median rule' in h:
            ITEMS.append(dict(rec, Source='Table S3', Model=S3rows[r[0]] + ', median activity rule'))
        else:
            s3_mean.append(dict(rec, Model=S3rows[r[0]]))
# Table S5: two fixed tie-break rules (hit counts, 1,000 permutations each, own nulls)
B = J(O + 'B_hits.json'); SN = {'HI': 'H. influenzae', 'SA': 'S. aureus', 'SPneu': 'S. pneumoniae', 'SPyo': 'S. pyogenes'}
for rule, lab in (('variance (as published)', 'largest variance, as published'), ('first index', 'first index')):
    for s, e in SN.items():
        p = B['tiebreak_fixed'][rule][s]['p']
        ITEMS.append(dict(Data='Monoterpenoids', Source='Table S5', Model='Route A1u hit count, tie-break ' + lab,
                          Endpoint=e, Metric='hit count', p=p, n=1000, k=k_exact(p, 1000)))
# Table S12: hit counts of the three routes (1,000 permutations)
for s, e in SN.items():
    for key, lab in (('A2', 'Route A2'), ('A1s', 'Route A1'), ('A1u', 'Route A1u')):
        p = B['calibration'][s][key]['p']
        ITEMS.append(dict(Data='Monoterpenoids', Source='Table S12', Model=lab + ' hit count', Endpoint=e, Metric='hit count',
                          p=p, n=1000, k=k_exact(p, 1000)))
# Table S13: random forest on the other pools (500) and PLS-DA on Morgan bits (2,000), exact values from AF2
hdr, rs = rows('Table1_new_all_arms'); ix = {h: i for i, h in enumerate(hdr)}
S13 = {'Random forest, out-of-bag, 182 representatives': 'Random forest, 182 label-free representatives',
       'Random forest, out-of-bag, Route B 12': 'Random forest, Route B 12 descriptors',
       'Morgan2 fingerprints, no selection': 'PLS-DA, Morgan2 bits, no selection',
       'Morgan2 fingerprints, nested filter': 'PLS-DA, Morgan2 bits, univariate filter nested'}
for r in rs:
    e = r[ix['Strain']]
    if e not in EP: continue
    lab = [v for pre, v in S13.items() if r[ix['Model']].startswith(pre)]
    if not lab: continue
    n = int(r[ix['Permutations']])
    for col, m in (('p_Q2', 'Q²'), ('p_AUC', 'AUC')):
        p = float(r[ix[col]])
        ITEMS.append(dict(Data='Monoterpenoids', Source='Table S13', Model=lab[0], Endpoint=e, Metric=m, p=p, n=n, k=k_exact(p, n)))
# random forest on Morgan bits (Table S13 a, last rows), 500 permutations
ITEMS += grid('Table_S13', 'Table S13', 'Monoterpenoids', 500,
              {'Morgan fingerprints (295 varying bits)': ('Random forest, Morgan2 bits', 'Q²')}, MONO)
hdr13, rs13 = rows('Table_S13')
for i, r in enumerate(rs13):
    if r[0].startswith('Morgan fingerprints (295'):
        nxt = rs13[i + 1]; assert nxt[0].startswith('AUC')
        for j, h in enumerate(hdr13[1:], 1):
            if h in EP:
                s = PAT.search(nxt[j]).group(1)
                ITEMS.append(dict(Data='Monoterpenoids', Source='Table S13', Model='Random forest, Morgan2 bits', Endpoint=h,
                                  Metric='AUC', p=float(s), n=500, k=k_rounded(s, 500)))
# Table S20: Gaucher, 2,000 permutations
hdr, rs = rows('Table_S20')
for r in rs:
    if r[0] in ('Nested Q² (p)', 'Nested AUC (p)'):
        m = 'Q²' if 'Q²' in r[0] else 'AUC'
        for j, h in enumerate(hdr[1:], 1):
            mm = PAT.search(r[j])
            if not mm: continue
            s = mm.group(1)
            ITEMS.append(dict(Data='Gaucher', Source='Table S20', Model='Route ' + h, Endpoint='Gaucher vs control', Metric=m,
                              p=float(s), n=2000, k=k_rounded(s, 2000)))
# new section: Routes A1 and A1u without the 3D descriptor blocks, 2,000 permutations
F = J(O + 'F_no3d_fixed.json')
for s, e in SN.items():
    for r in ('A1', 'A1u'):
        for key, m in (('p_q2', 'Q²'), ('p_auc', 'AUC')):
            p = F[s][r][key]
            ITEMS.append(dict(Data='Monoterpenoids', Source='3D section', Model='Route %s without 3D descriptors' % r,
                              Endpoint=e, Metric=m, p=p, n=2000, k=k_exact(p, 2000)))

# p-values reported only in the text of Additional file 1
CB = J(O + 'C_routeB.json')
for s, e in (('HI', 'H. influenzae'), ('SA', 'S. aureus'), ('SPneu', 'S. pneumoniae')):
    p = CB[s]['without_phenol_11']['p_q2']
    ITEMS.append(dict(Data='Monoterpenoids', Source='Route B section (text)', Model='Route B without the phenol indicator, 11 descriptors',
                      Endpoint=e, Metric='Q²', p=p, n=2000, k=k_exact(p, 2000)))
for s, e in (('SPneu', 'S. pneumoniae'), ('SPyo', 'S. pyogenes')):
    a = J(O + 'cache_topk_seeds/%s_A2_top15.json' % s); assert a['nperm'] == 2000
    ITEMS.append(dict(Data='Monoterpenoids', Source='Top-15 section (text)', Model='Route A2 top-15 by effect size, nested',
                      Endpoint=e, Metric='Q²', p=a['p_q2'], n=2000, k=k_exact(a['p_q2'], 2000)))

# ---------------------------------------------------------------- consistency checks against Tables 1 to 5
T1 = {(it['Model'], it['Endpoint'], it['Metric']): it['p'] for it in ITEMS if it['Source'] == 'Table 1'}
chk = []
for r in s3_mean:   # Table S3 mean-rule columns = Table 1
    p1 = T1[(r['Model'], r['Endpoint'], r['Metric'])]
    chk.append(('Table S3 mean rule', r['Model'], r['Endpoint'], r['Metric'], r['p'], p1, abs(r['p'] - p1) < 6e-5))
hdr, rs = rows('Table_S1')
for r in rs:
    for pre, mod in (('Q², 1 LV, grouping fixed', 'Route B, 12 descriptors, 1 LV'), ('Q², 3 LV, grouping fixed', 'Route B, 12 descriptors, 3 LV')):
        if r[0].startswith(pre):
            for j, h in enumerate(hdr[1:], 1):
                if h in EP:
                    pv = float(PAT.search(r[j]).group(1)); p1 = T1[(mod, h, 'Q²')]
                    chk.append(('Table S1 fixed', mod, h, 'Q²', pv, p1, abs(pv - p1) < 6e-5))
hdr, rs = rows('Table_S13')
for i, r in enumerate(rs):
    if r[0].startswith('All 2,208 descriptors'):
        for rr, m in ((r, 'Q²'), (rs[i + 1], 'AUC')):
            for j, h in enumerate(hdr[1:], 1):
                if h in EP:
                    pv = float(PAT.search(rr[j]).group(1)); p1 = T1[('Random forest, out-of-bag', h, m)]
                    chk.append(('Table S13 all descriptors', 'Random forest', h, m, pv, p1, abs(pv - p1) < 6e-5))
assert all(c[-1] for c in chk), [c for c in chk if not c[-1]]

# ---------------------------------------------------------------- intervals
seen = {}; dup = 0
for it in ITEMS:
    key = (it['Data'], it['Source'], it['Model'], it['Endpoint'], it['Metric'])
    assert key not in seen, key
    seen[key] = it
out = []
for it in ITEMS:
    lo, hi = cp(it['k'], it['n'])
    it.update(ci_low=lo, ci_high=hi, thresholds_in_ci=[t for t in THR if lo <= t <= hi])
    if it['thresholds_in_ci']: out.append(it)
by_src = {}
for it in ITEMS: by_src.setdefault(it['Source'], [0, 0]); by_src[it['Source']][0] += 1
for it in out: by_src[it['Source']][1] += 1
note = ('Table S5, 100 random tie-break rules: reported as median and range and as the share of rules with p <= 0.05 and '
        'p <= 0.0125. The per-rule p-values are not archived (B_hits.json keeps the summary only), so the share at 0.0167 '
        'requires re-running the counting step of jobB_hits.py with the same seeds.')
json.dump(dict(n_checked=len(ITEMS), n_in_table=len(out), thresholds=list(THR), by_source=by_src,
               identity_checks=[list(c) for c in chk], note=note, items=ITEMS, table_S10=out),
          open(OUT, 'w'), indent=1, ensure_ascii=False)
print('checked %d p-values, %d with 0.05 or 0.0167 in the 95 %% interval, %d identity checks passed' % (len(ITEMS), len(out), len(chk)))
for s, (a, b) in by_src.items(): print('  %-10s checked %3d  in Table S10 %2d' % (s, a, b))
print('--- entries outside Tables 1 to 5:')
for it in out:
    if it['Source'] in ('Table 1', 'Table 2', 'Table 3', 'Table 4', 'Table 5'): continue
    print('  %-10s %-62s %-14s %-9s p = %.4f  k = %3d/%d  CI %.4f-%.4f  (%s)' % (
        it['Source'], it['Model'], it['Endpoint'], it['Metric'], it['p'], it['k'], it['n'], it['ci_low'], it['ci_high'],
        ' and '.join('%.4g' % t for t in it['thresholds_in_ci'])))
print('NOTE:', note)
