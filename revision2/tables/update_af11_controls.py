# -*- coding: utf-8 -*-
"""Additional file 11: control arms of Table 1 and Table S14 added (response to Reviewer 2, point 6, text T6-1).
Eight rows per strain after the Route B row, from revision2/results/classification_controls.json, in the
order and with the labels of Table S19. The 35 existing rows are kept unchanged (checked).
Usage: python update_af11_controls.py IN.xlsx classification_controls.json OUT.xlsx"""
import sys, json, openpyxl
IN, CC, OUT = sys.argv[1:4]
R = json.load(open(CC))['results']
KEY = {'H. influenzae': 'HI', 'S. aureus': 'SA', 'S. pneumoniae': 'SPneu', 'S. pyogenes': 'SPyo', 'P. aeruginosa': 'PA'}
ARMS = [('No selection', 'PLS-DA, all 2,208 descriptors', 2208), ('No filter (182)', 'PLS-DA, 182 label-free representatives', 182),
        ('RF, all (OOB)', 'Random forest, all 2,208 descriptors', 2208), ('RF, 182 (OOB)', 'Random forest, 182 label-free representatives', 182),
        ('RF, B (OOB)', 'Random forest, Route B 12 descriptors', 12), ('RF, Morgan (OOB)', 'Random forest, Morgan bits', 295),
        ('Morgan, no selection', 'PLS-DA, Morgan bits, no selection', 295), ('Morgan, filter nested', 'PLS-DA, Morgan bits, filter nested', '–')]
wb = openpyxl.load_workbook(IN); ws = wb['Sheet1']
hdr = [c.value for c in ws[1]]
assert hdr == ['Strain', 'Model', 'Descriptors', 'TP', 'FN', 'FP', 'TN', 'Sensitivity', 'Specificity', 'BalancedAccuracy', 'MCC'], hdr
old = [list(r) for r in ws.iter_rows(min_row=2, values_only=True)]
assert len(old) == 35
new = []
for r in old:
    new.append(r)
    if str(r[1]).startswith('B (12'):
        k = KEY[r[0].replace(' †', '')]
        for lab, nm, nd in ARMS:
            v = R[k][nm]
            new.append([r[0], lab, nd, v['TP'], v['FN'], v['FP'], v['TN'], round(v['sens'], 3), round(v['spec'], 3), round(v['bal_acc'], 3), round(v['mcc'], 3)])
assert len(new) == 75
for i, row in enumerate(new):
    for j, v in enumerate(row):
        ws.cell(row=i + 2, column=j + 1, value=v)
info = wb['Info']
assert info['A1'].value.startswith('Additional file 11 - Classification metrics for every model reported in this study.')
info['A1'] = ('Additional file 11 - Classification metrics for every monoterpenoid model of Table 1 and Table S14: true and false positives '
              'and negatives, sensitivity, specificity, balanced accuracy and the Matthews correlation coefficient, at a decision threshold '
              'of 0.5 (a prediction of 0.5 or more counts as active).')
info['A2'] = 'Those of the second data set are given in Additional file 2 and those of the Gaucher data in Additional file 7.'
info['A3'] = 'Matches Table S19, whose legend explains the model labels.'
info['A4'] = None
wb.save(OUT)
chk = [list(r) for r in openpyxl.load_workbook(OUT)['Sheet1'].iter_rows(min_row=2, values_only=True)]
assert [r for r in chk if not str(r[1]).startswith(tuple(a[0] for a in ARMS))] == old, 'existing rows changed'
print('rows', len(chk), 'of which new', len(chk) - len(old))
