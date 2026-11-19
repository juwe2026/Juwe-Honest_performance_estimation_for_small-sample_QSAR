# -*- coding: utf-8 -*-
"""Additional file 2: changes that follow from the response to Reviewer 2 (approved 13.11.2026 16:31 and 16:44).
  Table_S4, Table_S11, Table_S19, Table_S20   rewritten cell by cell from Additional file 1 (T8-3, T8-4, T2-4/N2, T6/N3, T2-3)
  HitCounts_all_data_sets                     new columns: k, B and 95 % Clopper-Pearson interval of k/B for B = 1,000,
                                              and p, k and interval for B = 100,000 (revision2/analysis/hits_extended.py)
  Info                                        new row on these columns (approved 13.11.2026 18:11)
Usage: python update_af2_r2.py AF2_IN.xlsx AF1.docx HITS_DIR OUT.xlsx"""
import sys, json, copy, re, openpyxl, docx
from docx.oxml.ns import qn
IN, AF1, HX, OUT = sys.argv[1:5]
wb = openpyxl.load_workbook(IN)
d = docx.Document(AF1)
heads = {}; last = None; ti = 0
for el in d.element.body.iterchildren():
    if el.tag == qn('w:p'):
        m = re.match(r'(Table S\d+)[: ]', docx.text.paragraph.Paragraph(el, d).text)
        if m: last = m.group(1)
    elif el.tag == qn('w:tbl'):
        heads.setdefault(last, []).append(d.tables[ti]); ti += 1
def mirror(sheet, tab):
    ws = wb[sheet]; style_hdr = [copy.copy(c._style) for c in ws[1]]; style_body = [copy.copy(c._style) for c in ws[2]]
    ws.delete_rows(1, ws.max_row)
    for i, row in enumerate(tab.rows, 1):
        seen = []
        for c in row.cells:
            if c._tc in seen: continue
            seen.append(c._tc); j = len(seen)
            cell = ws.cell(i, j, c.text)
            st = (style_hdr if i == 1 else style_body)
            if j - 1 < len(st): cell._style = copy.copy(st[j - 1])
    return ws.max_row
LOG = []
for sheet, key, j in (('Table_S4', 'Table S4', 0), ('Table_S11', 'Table S11', 0), ('Table_S19', 'Table S19', 0), ('Table_S20', 'Table S20', 0)):
    n = mirror(sheet, heads[key][j]); LOG.append('%s: %d rows from Additional file 1' % (sheet, n))
# HitCounts_all_data_sets
mono = json.load(open(HX + '/hits_extended_monoterpenoids.json'))['results']
pep = json.load(open(HX + '/hits_extended_peptides.json'))['results']
gau = json.load(open(HX + '/hits_extended_gaucher.json'))['results']
KM = {'H. influenzae': 'H.I.', 'S. aureus': 'S.A.', 'S. pneumoniae': 'S.Pneu', 'S. pyogenes': 'S.Pyo', 'P. aeruginosa': 'P.A.'}
KP = {'E. coli ATCC 25922': 'EC', 'S. aureus ATCC 29213': 'SA', 'P. aeruginosa ATCC 27853': 'PA', 'A. baumannii ATCC 19606': 'AB'}
ws = wb['HitCounts_all_data_sets']; hdr = [c.value for c in ws[1]]
assert hdr[-1] == 'p' and 'p_B100000' not in hdr
NEW = ['k_B1000', 'CI95_low_k_over_B_B1000', 'CI95_high_k_over_B_B1000', 'p_B100000', 'k_B100000', 'CI95_low_k_over_B_B100000', 'CI95_high_k_over_B_B100000']
c0 = len(hdr)
for i, h in enumerate(NEW):
    c = ws.cell(1, c0 + 1 + i, h); c._style = copy.copy(ws.cell(1, c0)._style)
n = 0
for r in range(2, ws.max_row + 1):
    ds, ep, rt = ws.cell(r, 1).value, ws.cell(r, 2).value, ws.cell(r, 3).value.replace('Route ', '')
    if ds.startswith('Monoterpenoids'): v = mono[KM[ep]][rt]
    elif ds.startswith('Gaucher'): v = gau[rt]
    else: v = pep['A' if 'main' in ds else 'C'][KP[ep]][rt]
    assert abs(ws.cell(r, c0).value - v['B1000']['p']) < 6e-5, (ds, ep, rt)
    vals = [v['B1000']['k'], round(v['B1000']['ci95_k_over_B'][0], 5), round(v['B1000']['ci95_k_over_B'][1], 5),
            round(v['B100000']['p'], 6), v['B100000']['k'], round(v['B100000']['ci95_k_over_B'][0], 6), round(v['B100000']['ci95_k_over_B'][1], 6)]
    for i, x in enumerate(vals):
        c = ws.cell(r, c0 + 1 + i, x); c._style = copy.copy(ws.cell(r, c0)._style)
    n += 1
LOG.append('HitCounts_all_data_sets: %d rows, 7 new columns' % n)
info = wb['Info']; r = info.max_row + 1
info.cell(r, 1, 'HitCounts_all_data_sets, columns k_B1000 to CI95_high_k_over_B_B100000')._style = copy.copy(info.cell(r - 1, 1)._style)
info.cell(r, 2, 'k is the number of permuted counts at least as large as the observed one and B the number of permutations. '
             'The intervals are 95 % Clopper–Pearson intervals of k/B. The columns for B = 100,000 come from extending the same '
             'permutation streams, whose first 1,000 permutations are those of the column p (revision2/analysis/hits_extended.py).')._style = copy.copy(info.cell(r - 1, 2)._style)
LOG.append('Info: new row on the columns of HitCounts_all_data_sets')
wb.save(OUT)
print('\n'.join(LOG))
