# -*- coding: utf-8 -*-
"""results/revision/fallback_correction_tables.xlsx: column Script_entry_point next to the column Seed in the sheets
Fallback_correction and Fallback_correction_top15 (response to Reviewer 2, point 1). The script and entry point of
each seed are those of the sheet Seeds of Additional file 2. The sheet Info states this and the numbering of the
version of 6 October 2026. No other cell changes (checked).
Usage: python fallback_script_column.py IN.xlsx OUT.xlsx"""
import sys, copy, openpyxl
IN, OUT = sys.argv[1:3]
SCRIPT = {31: 'revision/jobI2_fallback.py / altorder.alt_nested_cv', 47: 'revision/jobI2_fallback.py / unsup_rep.unsup_nested_cv',
          7: 'revision/jobI2_fallback.py / nested_generic.nested_cv', 42: 'revision/jobL_a2hi.py / nested_generic.nested_cv',
          59: 'revision/jobJ2_topk.py / altorder.alt_nested_cv, 15 descriptors', 137: 'revision/jobJ2_topk.py / altorder.alt_viptop, one component',
          11: 'revision/jobJ2_topk.py / nested_generic.nested_cv_a2_topk, 15 descriptors', 23: 'revision/jobJ2_topk.py / a2vip.nested_cv_a2_viptop, two components'}
wb = openpyxl.load_workbook(IN)
before = {s: [[c.value for c in r] for r in wb[s].iter_rows()] for s in wb.sheetnames}
for s in ('Fallback_correction', 'Fallback_correction_top15'):
    ws = wb[s]; hdr = [c.value for c in ws[1]]; js = hdr.index('Seed') + 1; jn = ws.max_column + 1
    assert hdr[-1] == 'Seed'
    c = ws.cell(1, jn, 'Script_entry_point'); c._style = copy.copy(ws.cell(1, js)._style)
    for r in range(2, ws.max_row + 1):
        c = ws.cell(r, jn, SCRIPT[ws.cell(r, js).value]); c._style = copy.copy(ws.cell(r, js)._style)
info = wb['Info']; rows = {info.cell(r, 1).value: r for r in range(1, info.max_row + 1)}
assert info.cell(rows['Written by'], 2).value == 'revision2/tables/update_af2_rev2.py'
info.cell(rows['Written by'], 2).value = ('revision2/tables/update_af2_rev2.py. The column Script_entry_point was added by '
                                          'revision2/tables/fallback_script_column.py.')
r = info.max_row + 1; info.cell(r, 1, 'Numbering')._style = copy.copy(info.cell(r - 1, 1)._style)
info.cell(r, 2, 'Table and file numbers in this workbook are those of the version of 6 October 2026. Table 3 and Table 4 are now '
             'Tables S17 and S18, and Additional file 10 is now Additional file 15.')._style = copy.copy(info.cell(r - 1, 2)._style)
wb.save(OUT)
wb2 = openpyxl.load_workbook(OUT)
for s in wb2.sheetnames:
    after = [[c.value for c in r] for r in wb2[s].iter_rows()]
    if s == 'Info': continue
    for rb, ra in zip(before[s], after):
        assert ra[:len(rb)] == rb, s
print('column Script_entry_point added, no other cell changed')
