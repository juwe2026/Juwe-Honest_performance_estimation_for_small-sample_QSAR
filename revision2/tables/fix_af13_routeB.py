# -*- coding: utf-8 -*-
"""Additional file 13: Route B VIP columns joined by descriptor name (response to Reviewer 2, point 3).

In the sheets VIP_Lit_<strain> the columns RouteB_VIP_apparent and RouteB_VIP_nested_mean must hold
values only for the 12 descriptors of the Route B modelling set. The version of 6 November 2026
carried 9 further rows (18 cells) whose values were copies of other descriptors. This script empties
the Route B cells of every descriptor outside the modelling set, checks the remaining 12 values per
strain against figures/data_recalc/vip_routeB.json (key D12, the source of Figs. S2 and S3, rounded
to three decimals) and checks that no other cell of the workbook changes.
Usage: python fix_af13_routeB.py IN.xlsx OUT.xlsx
"""
import os as _osp, sys, json
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
import openpyxl

IN, OUT = sys.argv[1], sys.argv[2]
D12 = json.load(open(_osp.path.join(_ROOT, 'figures', 'data_recalc', 'vip_routeB.json')))['D12']
STRAINS = ['H.I.', 'S.A.', 'S.Pneu', 'S.Pyo', 'P.A.']
COLS = ('RouteB_VIP_apparent', 'RouteB_VIP_nested_mean')

wb = openpyxl.load_workbook(IN)
before = {s: [[c.value for c in r] for r in wb[s].iter_rows()] for s in wb.sheetnames}
log = []
for st in STRAINS:
    ws = wb['VIP_Lit_' + st]
    hdr = [c.value for c in ws[1]]
    ia, inn = hdr.index(COLS[0]), hdr.index(COLS[1])
    kept = {}
    for row in ws.iter_rows(min_row=2):
        name = row[0].value
        if not name:
            continue
        if name in D12[st]:
            a, n = round(D12[st][name]['apparent'], 3), round(D12[st][name]['nested'], 3)
            assert (row[ia].value, row[inn].value) == (a, n), (st, name, row[ia].value, row[inn].value, a, n)
            kept[name] = (a, n)
        elif row[ia].value is not None or row[inn].value is not None:
            log.append((st, name, row[ia].coordinate, row[ia].value, row[inn].coordinate, row[inn].value))
            row[ia].value = None; row[inn].value = None
    assert sorted(kept) == sorted(D12[st]) and len(kept) == 12, (st, sorted(kept))
wb.save(OUT)

wb2 = openpyxl.load_workbook(OUT)
changed = set((x[0], x[2]) for x in log) | set((x[0], x[4]) for x in log)
for s in wb2.sheetnames:
    after = [[c.value for c in r] for r in wb2[s].iter_rows()]
    for i, (rb, ra) in enumerate(zip(before[s], after)):
        for j, (vb, va) in enumerate(zip(rb, ra)):
            if vb != va:
                coord = openpyxl.utils.get_column_letter(j + 1) + str(i + 1)
                st = s.replace('VIP_Lit_', '')
                assert (st, coord) in changed, (s, coord, vb, va)
for st in STRAINS:
    ws = wb2['VIP_Lit_' + st]
    hdr = [c.value for c in ws[1]]
    n = sum(1 for r in ws.iter_rows(min_row=2, values_only=True) if r[0] and r[hdr.index(COLS[0])] is not None)
    assert n == 12, (st, n)
print('emptied %d rows (%d cells):' % (len(log), 2 * len(log)))
for x in log:
    print('  %-7s %-18s %s=%s %s=%s' % x)
print('every strain: exactly 12 Route B rows, all equal to vip_routeB.json D12, no other cell changed')
