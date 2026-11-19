# -*- coding: utf-8 -*-
"""Additional file 2: approved changes P2 (EW 06.11.2026 13:56 Berlin, 'P2 Option a').

Usage: python3 af2_p2.py AF2_IN.xlsx MAIN.docx AF2_OUT.xlsx "YYYY-MM-DD HH:MM CET"
P2a  new sheet Main_Table_3, written cell by cell from Table 3 of the main text as printed (minus signs as ASCII
     hyphens, as in the other table sheets), formatted like Main_Table_2, placed after it. Info row 'Table sheets'.
P2b  Info row 'Number_check'.
P2c  Info row 'Generated': the time of writing."""
import copy, re, sys
import openpyxl, docx
from docx.oxml.ns import qn

IN, MAIN, OUT, WHEN = sys.argv[1:5]
assert re.fullmatch(r'\d{4}-\d\d-\d\d \d\d:\d\d CET', WHEN)
wb = openpyxl.load_workbook(IN)
assert 'Main_Table_3' not in wb.sheetnames

# main text Table 3: the table that follows the caption 'Table 3'
d = docx.Document(MAIN); last = None; ti = 0; tab = None
for el in d.element.body.iterchildren():
    if el.tag == qn('w:p'):
        m = re.match(r'(Table \d+)\s', ''.join(t.text or '' for t in el.iter(qn('w:t'))))
        if m: last = m.group(1)
    elif el.tag == qn('w:tbl'):
        if last == 'Table 3': tab = d.tables[ti]
        ti += 1
assert tab is not None
rows = []
for r in tab.rows:
    seen, vals = [], []
    for c in r.cells:
        if c._tc in seen: continue
        seen.append(c._tc); vals.append(c.text.replace('\u2212', '-'))
    rows.append(vals)
assert len(rows) == 5 and all(len(r) == 8 for r in rows) and rows[0][0] == 'Hypothesis'

tpl = wb['Main_Table_2']
ws = wb.create_sheet('Main_Table_3', index=wb.sheetnames.index('Main_Table_2') + 1)
for i, r in enumerate(rows, 1):
    for j, v in enumerate(r, 1):
        c = ws.cell(i, j, v if v != '' else None)
        src = tpl.cell(1 if i == 1 else 2, 2 if j > 1 else 1)
        c.font = copy.copy(src.font); c.alignment = copy.copy(src.alignment); c.border = copy.copy(src.border)
        c.fill = copy.copy(src.fill); c.number_format = src.number_format
for k, dim in tpl.column_dimensions.items():
    ws.column_dimensions[k].width = dim.width

info = wb['Info']
def setrow(key, old, new, whole=False):
    hit = [r for r in range(1, info.max_row + 1) if info.cell(r, 1).value == key]
    assert len(hit) == 1, key
    c = info.cell(hit[0], 2)
    if whole:
        c.value = new; return
    assert c.value.count(old) == 1, (key, c.value[:80])
    c.value = c.value.replace(old, new)
setrow('Table sheets',
       'Main_Table_1 to Main_Table_5 and Table_S1 to Table_S22 reproduce every table of the main text and of Additional file 1 cell by cell, as printed, except that Table_S11 gives',
       'Main_Table_1 to Main_Table_3 and Table_S1 to Table_S26 reproduce every table of the main text and of Additional file 1 cell by cell, as printed, except that Table_S13 gives')
setrow('Number_check', 'every value and p-value printed in Tables 1 to 5, checked against',
       'every value and p-value printed in Tables 1 and 2 of the main text and in Tables S16 to S18 of Additional file 1, checked against')
old_gen = [info.cell(r, 2).value for r in range(1, info.max_row + 1) if info.cell(r, 1).value == 'Generated'][0]
assert old_gen == '2026-10-06 11:29 CEST', old_gen
setrow('Generated', None, WHEN, whole=True)
wb.save(OUT)
print('Main_Table_3', len(rows), 'rows x', len(rows[0]), '| Info rows Table sheets, Number_check, Generated set |', wb.sheetnames[-4:])
