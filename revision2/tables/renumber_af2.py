# -*- coding: utf-8 -*-
"""Additional file 2 in the approved numbering (EW 06.11.2026 13:34 Berlin, N2).

Usage: python3 renumber_af2.py AF2_IN.xlsx AF2_OUT.xlsx

1. Sheets Table_S<n>(_part2) and Main_Table_5 renamed to the new numbers, table sheets ordered by number.
2. References in all text cells renumbered (nummern_freigegeben.rn_old_text: former main Tables 2 to 4 become
   Tables S16 to S18). The Info rows 'Table sheets' and 'Number_check' are left for their approved wording.
Values that are not text are not touched (checked)."""
import os, re, sys, json
import openpyxl
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nummern as N
import nummern_freigegeben as NF

IN, OUT = sys.argv[1:3]
wb = openpyxl.load_workbook(IN)
before = {ws.title: [[c.value for c in row] for row in ws.iter_rows()] for ws in wb.worksheets}
SKIP = {('Info', 'Table sheets'), ('Info', 'Number_check')}
changed = []
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if not isinstance(c.value, str): continue
            if ws.title == 'Info' and c.column == 2 and (ws.title, ws.cell(c.row, 1).value) in SKIP: continue
            v = NF.rn_old_text(c.value)
            if v != c.value: changed.append((ws.title, c.coordinate, c.value[:70], v[:70])); c.value = v
# rename sheets
ren = {}
for name in wb.sheetnames:
    m = re.match(r'Table_S(\d+)(_part2)?$', name)
    if m: ren[name] = 'Table_S%d%s' % (NF.TS[int(m.group(1))], m.group(2) or '')
    elif name == 'Main_Table_5': ren[name] = 'Main_Table_2'
for a in ren: wb[a].title = 'tmp_' + a
for a, b in ren.items(): wb['tmp_' + a].title = b
key = lambda n: (0, 0, '') if not (n.startswith('Table_S') or n.startswith('Main_Table_')) else \
      ((1, int(re.findall(r'\d+', n)[0]), n) if n.startswith('Table_S') else (2, int(n.split('_')[-1]), n))
order = sorted(range(len(wb._sheets)), key=lambda i: (key(wb._sheets[i].title)[0] if key(wb._sheets[i].title)[0] else 0,
                                                       i if key(wb._sheets[i].title)[0] == 0 else 0, key(wb._sheets[i].title)[1:]))
wb._sheets = [wb._sheets[i] for i in order]
wb.save(OUT)
# non-text cells unchanged
wb2 = openpyxl.load_workbook(OUT)
inv = {b: a for a, b in ren.items()}
nontext = 0
for ws in wb2.worksheets:
    old = before[inv.get(ws.title, ws.title)]
    new = [[c.value for c in row] for row in ws.iter_rows()]
    assert len(old) == len(new)
    for ro, rw in zip(old, new):
        for a, b in zip(ro, rw):
            if not isinstance(a, str): assert a == b or (a != a and b != b), (ws.title, a, b); nontext += 1
json.dump({'renamed': ren, 'cells_changed': changed, 'sheet_order': wb2.sheetnames}, open(os.path.splitext(OUT)[0] + '_log.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
print('sheets renamed', len(ren), '| text cells changed', len(changed), '| non-text cells identical', nontext)
print(wb2.sheetnames)
