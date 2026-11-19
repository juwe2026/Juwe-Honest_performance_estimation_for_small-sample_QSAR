# -*- coding: utf-8 -*-
"""Additional file 2 (Revision 2, working numbers): content updates that follow from approved changes.

Inputs: af2_rev1.xlsx (Additional file 2 as submitted), S10.json (cp_borderline_all_rev2.py, 325 values),
TieBreak.json (jobB_tiebreak_0167.py, per-rule p-values), AF1 docx of Revision 2 (for Table S23 cell by cell).
Changes:
  Table_S5         last row: p <= 0.0167 shares, derived from TieBreak.json (approved 06.11.2026)
  Table_S8_part2   semicolons in cells replaced by commas (approved 06.11.2026)
  Table_S10        rebuilt: 67 rows, thresholds 0.05 and 0.0167 (approved 06.11.2026)
  Table_S19        new row 'Four evaluable endpoints (without P. aeruginosa)' (approved 05.11.2026)
  Endpoint_correlation  new block for the four evaluable monoterpenoid endpoints, computed from Activity_vectors_31 (approved 05.11.2026)
  MonteCarlo_uncertainty  data rows replaced by the 67 rows of Table S10, unrounded (illustration rows kept)
  MonteCarlo_all_checked  NEW: all 325 checked p-values with k, n and interval
  TieBreak_sensitivity  summary rows of the 100 random rules state the share at 0.0167
  TieBreak_random_rules  NEW: the 100 random tie-break rules, one row per rule and endpoint
  Table_S23        NEW: copied cell by cell from Additional file 1
  Main_Table_2/3/4 renamed Table_S24/S25/S26 (contents unchanged)
  approved 06.11.2026 08:28: D1 fallback sheets to the archive (seed column corrected), D2 Seeds, D3 semicolons and Info,
  D4 English column heads of Null_distribution_audit, D5 provenance sentence, D6 peptide eigenvalues, E10 labels
  approved 06.11.2026 08:43: D7 n_tests of the Gaucher rows, D8 Seeds F32, D9 wording of the Table S5 row label
Usage: python update_af2_rev2.py AF1_REV2.docx OUT.xlsx ARCHIVE_DIR FALLBACK_OUT.xlsx
  (ARCHIVE_DIR: unpacked code archive with results/revision/, aenderungen_af2.py next to this script)
"""
import json, sys, copy, os, glob
import openpyxl, docx
from openpyxl.styles import Font

AF1, OUT, ARCH, FB_OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aenderungen_af2 as AE                         # approved text changes D2, D3, D4, E10 (06.11.2026 08:28)
wb = openpyxl.load_workbook('af2_rev1.xlsx')
LOG = []
BOLD = Font(name='Calibri', size=11, bold=True)

def header(ws, cols):
    for k, c in enumerate(cols, 1):
        ws.cell(row=1, column=k, value=c).font = BOLD

# Table_S5
ws = wb['Table_S5']
assert ws['A6'].value == 'share of rules with p ≤ 0.05 / p ≤ 0.0125'
assert [ws.cell(6, k).value for k in range(2, 7)] == ['1 % / 0 %', '100 % / 100 %', '100 % / 57 %', '0 % / 0 %', '100 % / 90 %']
TB = json.load(open('TieBreak.json')); assert TB['reproduces_B_hits']
ws['A6'] = 'share of rules with p ≤ 0.05 / p ≤ 0.0167'
for col, key in zip('BCDEF', ('HI', 'SA', 'SPneu', 'SPyo', 'PA')):      # derived from TieBreak.json, not typed
    s = TB['summary'][key]; ws[col + '6'] = '%.0f %% / %.0f %%' % (100 * s['frac_p_le_05'], 100 * s['frac_p_le_0167'])
assert [ws.cell(6, k).value for k in range(2, 7)] == ['1 % / 0 %', '100 % / 100 %', '100 % / 83 %', '0 % / 0 %', '100 % / 100 %']
LOG.append('Table_S5: last row (from TieBreak.json)')
# Table_S8_part2
ws = wb['Table_S8_part2']; n = 0
for row in ws.iter_rows():
    for c in row:
        if isinstance(c.value, str) and ';' in c.value:
            assert c.value.count('; ') == c.value.count(';'); c.value = c.value.replace('; ', ', '); n += 1
LOG.append('Table_S8_part2: %d cells semicolon -> comma' % n)
# Table_S10
S = json.load(open('S10.json')); rows = S['table_S10']
order = ['Table 1', 'Table 2', 'Table 3', 'Table 4', 'Table 5', 'Table S1', 'Table S2', 'Table S3', 'Table S5', 'Table S12',
         'Table S13', 'Table S20', '3D section', 'Route B section (text)', 'Top-15 section (text)']
lab = {'Table 2': 'Table S24', 'Table 3': 'Table S25', 'Table 4': 'Table S26', '3D section': 'Table S23',
       'Route B section (text)': 'Route B in detail (text)', 'Top-15 section (text)': 'Fixed model size (text)'}
rows.sort(key=lambda r: order.index(r['Source']))
ws = wb['Table_S10']; hdr = [ws.cell(1, k).value for k in range(1, 9)]
assert hdr == ['Source', 'Model', 'Endpoint', 'Metric', 'p', 'k / n', '95 % CI of k/n', 'Threshold in CI']
ws.delete_rows(2, ws.max_row)
for i, r in enumerate(rows, 2):
    ep = r['Endpoint'] + (' (peptides)' if r['Data'] == 'Octapeptides' else '')
    vals = [lab.get(r['Source'], r['Source']), r['Model'].replace(*AE.LABEL_S5), ep, r['Metric'], '%.4f' % ((r['k'] + 1) / (r['n'] + 1)),
            '%d / %d' % (r['k'], r['n']), '%.4f–%.4f' % (r['ci_low'], r['ci_high']),
            ' and '.join('0.05' if t == 0.05 else '0.0167' for t in r['thresholds_in_ci'])]
    for k, v in enumerate(vals, 1): ws.cell(i, k, v)
LOG.append('Table_S10: %d rows' % len(rows))
# MonteCarlo_uncertainty: keep the illustration rows, replace the data rows
ws = wb['MonteCarlo_uncertainty']
assert [ws.cell(1, k).value for k in range(1, 8)] == ['Kind', 'Permutations', 'k', 'p_reported', 'CI95_low', 'CI95_high', 'Smallest_attainable_p']
last_ill = max(r for r in range(2, ws.max_row + 1) if ws.cell(r, 1).value == 'illustration')
ws.delete_rows(last_ill + 1, ws.max_row)
for i, r in enumerate(rows, last_ill + 1):
    kind = '%s, %s, %s%s, %s' % (lab.get(r['Source'], r['Source']), r['Model'].replace(*AE.LABEL_S5), r['Endpoint'], ' (peptides)' if r['Data'] == 'Octapeptides' else '', r['Metric'])
    vals = [kind, r['n'], r['k'], round((r['k'] + 1) / (r['n'] + 1), 5), round(r['ci_low'], 5), round(r['ci_high'], 5), round(1 / (r['n'] + 1), 5)]
    for k, v in enumerate(vals, 1): ws.cell(i, k, v)
LOG.append('MonteCarlo_uncertainty: %d data rows' % len(rows))
# MonteCarlo_all_checked (new)
ws = wb.create_sheet('MonteCarlo_all_checked', index=wb.sheetnames.index('MonteCarlo_uncertainty') + 1)
header(ws, ['Data_set', 'Source', 'Model', 'Endpoint', 'Metric', 'Permutations', 'k', 'p_reported', 'CI95_low', 'CI95_high', 'Interval_contains'])
for i, r in enumerate(sorted(S['items'], key=lambda r: order.index(r['Source'])), 2):
    vals = [r['Data'], lab.get(r['Source'], r['Source']), r['Model'].replace(*AE.LABEL_S5), r['Endpoint'], r['Metric'], r['n'], r['k'],
            round((r['k'] + 1) / (r['n'] + 1), 5), round(r['ci_low'], 5), round(r['ci_high'], 5),
            ' and '.join('0.05' if t == 0.05 else '0.0167' for t in r['thresholds_in_ci'])]
    for k, v in enumerate(vals, 1): ws.cell(i, k, v)
assert ws.max_row - 1 == S['n_checked'] == 325
LOG.append('MonteCarlo_all_checked: %d rows (new)' % (ws.max_row - 1))
# TieBreak_sensitivity and TieBreak_random_rules
NM = {'HI': 'H. influenzae', 'SA': 'S. aureus', 'SPneu': 'S. pneumoniae', 'SPyo': 'S. pyogenes', 'PA': 'P. aeruginosa'}
ws = wb['TieBreak_sensitivity']
for r in range(2, ws.max_row + 1):
    if ws.cell(r, 1).value == '100 random tie-break rules':
        key = [k for k, v in NM.items() if v == ws.cell(r, 3).value][0]; s = TB['summary'][key]
        old = ws.cell(r, 7).value
        new = '%.4f-%.4f (median %.4f), p<=0.05 in %.0f%% of rules, p<=0.0167 in %.0f%%' % (s['p_min'], s['p_max'], s['p_median'], 100 * s['frac_p_le_05'], 100 * s['frac_p_le_0167'])
        assert old.split(';')[0] == new.split(',')[0], (old, new)
        ws.cell(r, 7, new)
ws = wb.create_sheet('TieBreak_random_rules', index=wb.sheetnames.index('TieBreak_sensitivity') + 1)
header(ws, ['Rule', 'n_groups', 'Endpoint', 'Observed_hits', 'Null_mean', 'Null_SD', 'p'])
i = 2
for rule in sorted(TB['per_rule']):
    pr = TB['per_rule'][rule]
    for key in ('HI', 'SA', 'SPneu', 'SPyo', 'PA'):
        x = pr[key]
        for k, v in enumerate([rule, pr['n_groups'], NM[key], x['obs'], round(x['null_mean'], 3), round(x['null_sd'], 3), round(x['p'], 4)], 1): ws.cell(i, k, v)
        i += 1
LOG.append('TieBreak: summary 0.0167, TieBreak_random_rules %d rows (new)' % (i - 2))
# Table_S19 and Endpoint_correlation
ws = wb['Table_S19']
r8 = [r for r in range(1, ws.max_row + 1) if ws.cell(r, 1).value == 'Effective number of endpoints, Li and Ji / Nyholt']
assert len(r8) == 1 and ws.cell(r8[0], 2).value == '3.0 / 3.4'
# four evaluable endpoints: phi, eigenvalues and effective numbers computed from the sheet Activity_vectors_31 (not typed)
import numpy as np
av = wb['Activity_vectors_31']
names = ['H. influenzae', 'S. aureus', 'S. pneumoniae', 'S. pyogenes']
avc = {av.cell(1, k).value: k for k in range(1, av.max_column + 1)}
assert av.max_row - 1 == 31
Lm = np.array([[float(av.cell(r, avc[nm]).value) for r in range(2, av.max_row + 1)] for nm in names])
phi = np.corrcoef(Lm); ev = np.sort(np.linalg.eigvalsh(phi))[::-1]; m = len(ev)
li_ji = float(sum((e >= 1) + (e - np.floor(e)) for e in np.abs(ev)))          # Li and Ji 2005
nyholt = float(1 + (m - 1) * (1 - np.var(ev, ddof=1) / m))                     # Nyholt 2004, Cheverud 2001
L5 = np.array([[float(av.cell(r, avc[nm]).value) for r in range(2, av.max_row + 1)] for nm in names + ['P. aeruginosa']])
p5 = np.corrcoef(L5)                                                          # cross-check against the existing five-endpoint block
ec = wb['Endpoint_correlation']
for j in range(5):
    for k in range(5): assert abs(ec.cell(2 + j, 3 + k).value - round(p5[j, k], 4)) < 1e-9
ws.insert_rows(r8[0] + 1)
ws.cell(r8[0] + 1, 1, 'Four evaluable endpoints (without P. aeruginosa), Li and Ji / Nyholt'); ws.cell(r8[0] + 1, 2, '%.1f / %.1f' % (li_ji, nyholt))
assert ws.cell(r8[0] + 1, 2).value == '3.0 / 2.7'                             # as approved for Table S19
ws = wb['Endpoint_correlation']
cols = {ws.cell(1, k).value: k for k in range(1, ws.max_column + 1)}
blk = 'four evaluable monoterpenoid endpoints (without P. aeruginosa)'
ins = 9                                              # after the block of the five monoterpenoid endpoints (rows 2 to 8)
assert ws.cell(8, 2).value.startswith('effective number of endpoints') and ws.cell(9, 1).value.startswith('37 peptides')
ws.insert_rows(ins, amount=6)
for j, nm in enumerate(names):
    ws.cell(ins + j, 1, '31 monoterpenoids, ' + blk); ws.cell(ins + j, 2, nm)
    for m, nm2 in enumerate(names): ws.cell(ins + j, cols['phi vs ' + nm2], round(float(phi[j, m]), 4))
ws.cell(ins + 4, 1, '31 monoterpenoids, ' + blk); ws.cell(ins + 4, 2, 'eigenvalues (descending)')
for m, v in enumerate(ev): ws.cell(ins + 4, cols['phi vs ' + names[m]], round(float(v), 4))
ws.cell(ins + 5, 1, '31 monoterpenoids, ' + blk); ws.cell(ins + 5, 2, 'effective number of endpoints: Li and Ji %.2f, Nyholt/Cheverud %.2f' % (li_ji, nyholt))
LOG.append('Table_S19: new row, Endpoint_correlation: new block (approved 05.11.2026)')
# Table_S23 from Additional file 1 (cell by cell)
d = docx.Document(AF1)
t23 = [t for t in d.tables if t.rows[1].cells[0].text.startswith('(a) Nested Q²') and 'without 3D' in t.rows[3].cells[0].text]
assert len(t23) == 1
ws = wb.create_sheet('Table_S23', index=wb.sheetnames.index('Table_S22') + 1)
for i, row in enumerate(t23[0].rows, 1):
    for k, c in enumerate(row.cells, 1): ws.cell(i, k, c.text if (i > 1 or k > 1) else 'col0')
for k in range(1, 7): ws.cell(1, k).font = BOLD
LOG.append('Table_S23: %d rows from Additional file 1 (new)' % len(t23[0].rows))
# Main_Table_2/3/4 -> Table_S24/S25/S26, placed after Table_S23
for old, new in (('Main_Table_2', 'Table_S24'), ('Main_Table_3', 'Table_S25'), ('Main_Table_4', 'Table_S26')):
    wb[old].title = new
    wb.move_sheet(new, offset=wb.sheetnames.index('Table_S23') + (int(new[-2:]) - 23) - wb.sheetnames.index(new))
LOG.append('Main_Table_2/3/4 -> Table_S24/S25/S26')

# ================================================================ approved 06.11.2026 08:28: D1 to D6 and E10
def setc(ws, coord, old, new):
    v = ws[coord].value
    if old is None: assert str(v).startswith('archived file'), (ws.title, coord, v)
    else: assert v == old, (ws.title, coord, v)
    ws[coord].value = new
# D1: the four sheets on the fallback correction go to the code archive, with the seed column corrected
FB = ['Fallback_correction', 'Fallback_correction_top15', 'Fallback_correction_TableS2', 'Fallback_affected_files']
fb = openpyxl.load_workbook('af2_rev1.xlsx')
for nm in fb.sheetnames:
    if nm not in FB: del fb[nm]
RV = os.path.join(ARCH, 'results', 'revision')
KEY = {'H. influenzae': 'HI', 'S. aureus': 'SA', 'S. pneumoniae': 'SPneu', 'S. pyogenes': 'SPyo', 'P. aeruginosa': 'PA'}
IF = json.load(open(os.path.join(RV, 'I_fallback_seeds.json'))); HI2 = json.load(open(os.path.join(RV, 'A2_HI_seed42.json')))
ws = fb['Fallback_correction']; hdr = [c.value for c in ws[1]]; assert hdr[10] == 'Seed'
SEED_T1 = {'A1': 31, 'A1u': 47, 'A2': 7}
for r in range(2, ws.max_row + 1):
    ep, route = ws.cell(r, 1).value, ws.cell(r, 2).value
    src = HI2 if (route == 'A2' and ep == 'H. influenzae') else IF[KEY[ep]][route]
    assert round(src['p_q2'], 4) == ws.cell(r, 5).value and round(src['p_auc'], 4) == ws.cell(r, 8).value, (ep, route)
    assert ws.cell(r, 11).value == 702
    ws.cell(r, 11).value = 42 if (route == 'A2' and ep == 'H. influenzae') else SEED_T1[route]
ws = fb['Fallback_correction_top15']; hdr = [c.value for c in ws[1]]; assert hdr[15] == 'Seed'
ARM = [('Route A1, top 15 by effect', 'A1_top15', 59), ('Route A1, top 15 by VIP', 'A1_viptop1', 137),
       ('Route A2, top 15 by effect', 'A2_top15', 11), ('Route A2, top 15 by VIP', 'A2_viptop', 23)]
for r in range(2, ws.max_row + 1):
    arm = [a for a in ARM if ws.cell(r, 2).value.startswith(a[0])]; assert len(arm) == 1
    src = json.load(open(os.path.join(RV, 'cache_topk_seeds', '%s_%s.json' % (KEY[ws.cell(r, 3).value], arm[0][1]))))
    assert (round(src['q2'], 4), round(src['auc'], 4), round(src['p_q2'], 4), round(src['p_auc'], 4), src['nperm']) == \
           (ws.cell(r, 5).value, ws.cell(r, 7).value, ws.cell(r, 9).value, ws.cell(r, 11).value, ws.cell(r, 13).value), r
    assert ws.cell(r, 16).value == 702
    ws.cell(r, 16).value = arm[0][2]
inf = fb.create_sheet('Info', 0)
for i, (a, b) in enumerate([('Field', 'Value'),
        ('Content', 'the four sheets on the correction of the empty-selection fallback that were part of Additional file 2 of the version submitted on 06.10.2026, with the values before and after the correction'),
        ('Seed column', 'gives the seed of the run that produced the values. In the submitted Additional file 2 this column showed 702 in the sheets Fallback_correction and Fallback_correction_top15, the seed of a first run that was discarded. The values themselves were from the runs with the seeds given here, checked against I_fallback_seeds.json, A2_HI_seed42.json and cache_topk_seeds/'),
        ('Written by', 'revision2/tables/update_af2_rev2.py')], 1):
    inf.cell(i, 1, a); inf.cell(i, 2, b)
    if i == 1: inf.cell(i, 1).font = BOLD; inf.cell(i, 2).font = BOLD
fb.save(FB_OUT)
for nm in FB: del wb[nm]
LOG.append('D1: 4 sheets moved to the archive (%s), seed column 35 rows corrected and checked against the result objects' % os.path.basename(FB_OUT))
# D2: Seeds
ws = wb['Seeds']
for sh, c, old, new in AE.SEEDS: setc(ws, c, old, new)
assert ws.max_row == 32
for k, v in enumerate(AE.SEEDS_NEW_ROW, 1):
    ws.cell(33, k, v); ws.cell(33, k).font = copy.copy(ws.cell(32, k).font); ws.cell(33, k).alignment = copy.copy(ws.cell(32, k).alignment)
LOG.append('D2: seeds %d cells, new row 33' % len(AE.SEEDS))
# D3: Info and repeated cells, D4: column heads of Null_distribution_audit
for sh, c, old, new in AE.INFO + AE.NULLAUDIT: setc(wb[sh], c, old, new)
for sh, col, old, new in AE.REPEATED:
    ws = wb[sh]; n = 0
    for r in range(1, ws.max_row + 1):
        if ws['%s%d' % (col, r)].value == old: ws['%s%d' % (col, r)].value = new; n += 1
    assert n > 0, (sh, old)
    LOG.append('D3/E10: %s column %s, %d cells' % (sh, col, n))
for old, new in AE.SHEET_RENAME: wb[old].title = new
# E10: row head of Table_S5
setc(wb['Table_S5'], 'A2', AE.TABLE_S5_A2[0], AE.TABLE_S5_A2[1])
# D5: provenance sentence (check_af2.py compares every sheet with its source)
setc(wb['Info'], 'B5', 'the result sheets are written directly from the deposited result objects, and the table sheets from the manuscript files as printed (minus signs as ASCII hyphens), both by mkworkbook.py; no value is typed by hand',
     'the result sheets are taken by script from the result objects deposited in Additional file 7, and the table sheets from the main text and Additional file 1 as printed (minus signs as ASCII hyphens). No value is typed by hand. check_af2.py in Additional file 7 compares every sheet with its source')
# D7 (approved 08:43): number of tests for the Gaucher rows of HitCounts_all_data_sets, from G_gaucher_hits.json
ws = wb['HitCounts_all_data_sets']; hdr = [c.value for c in ws[1]]; GH = json.load(open(os.path.join(RV, 'G_gaucher_hits.json')))
cn, cr, cs, cnom = hdr.index('n_tests') + 1, hdr.index('Route') + 1, hdr.index('Null_SD') + 1, hdr.index('Inflation_vs_nominal_5pct') + 1
n7 = 0
for r in range(2, ws.max_row + 1):
    if str(ws.cell(r, 1).value).startswith('Gaucher'):
        rt = ws.cell(r, cr).value.replace('Route ', ''); nt = GH['n_var'] if rt == 'A2' else GH['n_groups']
        assert ws.cell(r, cn).value is None
        # the printed inflation factor is null SD / binomial SD at 5 %, which holds only with this number of tests
        assert abs(ws.cell(r, cs).value / (nt * 0.05 * 0.95) ** 0.5 - ws.cell(r, cnom).value) < 0.001, (rt, nt)
        ws.cell(r, cn).value = nt; n7 += 1
assert n7 == 3
LOG.append('D7: n_tests of the 3 Gaucher rows (590, 428, 428)')
# D6: eigenvalues of the peptide blocks so that they descend from left to right
ws = wb['Endpoint_correlation']; EC = json.load(open(os.path.join(RV, 'endpoint_corr.json')))
for r in range(2, ws.max_row + 1):
    a, b = ws.cell(r, 1).value, ws.cell(r, 2).value
    if a.startswith('37 peptides') and b == 'eigenvalues (descending)':
        cols = [k for k in range(3, ws.max_column + 1) if ws.cell(r, k).value is not None]
        vals = sorted((ws.cell(r, k).value for k in cols), reverse=True)
        ref = EC['lone_rule_A' if 'rule A' in a else 'lone_rule_C']['eigen']
        assert vals == [round(float(x), 4) for x in ref], (a, vals, ref)
        for k, v in zip(cols, vals): ws.cell(r, k).value = v
        LOG.append('D6: %s eigenvalues descending' % a)
wb.save(OUT)
print('\n'.join(LOG)); print(wb.sheetnames)
