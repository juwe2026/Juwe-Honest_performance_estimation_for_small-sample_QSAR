# -*- coding: utf-8 -*-
"""Applies aenderungen_af3_13 to the submitted Additional files 3, 4, 6 and 8 to 13 and checks the result.

Usage: python3 apply_af3_13.py SRC_DIR OUT_DIR TIMESTAMP
Writes the renumbered files 'Additional_file_<new>_<name> TIMESTAMP.xlsx', a change list (JSON) and prints the checks:
every number cell unchanged, no semicolon left outside SMARTS patterns, no reference to an earlier version left,
every reference to a table, figure or additional file valid in the new numbering."""
import json, os, re, sys, copy
import openpyxl
from openpyxl.utils import column_index_from_string
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import aenderungen_af3_13 as A
import nummern as N
import nummern_freigegeben as NF

SRC, OUT, TS = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
HIST = re.compile(r'submitted|earlier|previous|revision|revised|REVISION|no longer|recomputed|regenerated|corrected value|the deposited code|original computation', re.I)
SMARTS = lambda v: bool(re.search(r'\[[#A-Za-z0-9;$!()=@+-]+\]', v)) and ';' in re.sub(r'\[[^\]]*\]', '', v) is False
changes, report = [], []
for old, F in A.FILES.items():
    wb = openpyxl.load_workbook(os.path.join(SRC, F['src']))
    before = {ws.title: {c.coordinate: c.value for row in ws.iter_rows() for c in row if c.value is not None} for ws in wb.worksheets}
    done = set()
    for (o, sh, col, head) in A.DELETE_COLUMNS:
        if o != old: continue
        ws = wb[sh]; k = column_index_from_string(col)
        assert ws.cell(1, k).value == head and k == ws.max_column, (sh, col)
        changes.append((old, sh, col + ':' + col, 'H', head + ' (column, %d values)' % (ws.max_row - 1), None))
        ws.delete_cols(k)
    for (o, sh, cell), (kind, new) in A.CELLS.items():
        if o != old: continue
        c = wb[sh][cell]; assert isinstance(c.value, str), (old, sh, cell)
        changes.append((old, sh, cell, kind, c.value, new)); c.value = new; done.add((sh, cell))
    for (o, prefix, cells) in A.SEMI_TO_COMMA:
        if o != old: continue
        for ws in wb.worksheets:
            if not ws.title.startswith(prefix): continue
            for cell in cells:
                v = ws[cell].value; assert v.count('; ') == 1, (ws.title, cell, v)
                nv = NF.rn_old_text(v.replace('; ', ', '))
                changes.append((old, ws.title, cell, 'S', v, nv)); ws[cell].value = nv; done.add((ws.title, cell))
    for (o, sh), (r1, r2, kind, lines) in A.BLOCKS.items():
        if o != old: continue
        ws = wb[sh]; assert ws.max_row == r2, (old, sh, ws.max_row)
        oldtxt = [ws.cell(r, 1).value for r in range(r1, r2 + 1)]
        assert all(ws.cell(r, k).value is None for r in range(r1, r2 + 1) for k in range(2, ws.max_column + 1))
        changes.append((old, sh, 'A%d:A%d' % (r1, r2), kind, '\n'.join(x or '' for x in oldtxt), '\n'.join(lines)))
        for i, line in enumerate(lines): ws.cell(r1 + i, 1).value = line if line else None
        n_extra = (r2 - r1 + 1) - len(lines)
        assert n_extra >= 0
        if n_extra: ws.delete_rows(r1 + len(lines), n_extra)
        for i in range(len(lines)): done.add((sh, 'A%d' % (r1 + i)))
    ren = A.SHEET_RENAMES.get(old, {})
    for ws in list(wb.worksheets):
        for a, b in ren.items():
            if ws.title == a or (a.endswith('_') and ws.title.startswith(a)):
                nt = b + ws.title[len(a):] if a.endswith('_') else b
                changes.append((old, ws.title, '(sheet name)', 'M', ws.title, nt)); ws.title = nt
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row:
                if not isinstance(c.value, str) or (ws.title, c.coordinate) in done: continue
                v = c.value
                for a, b in ren.items(): v = v.replace(a, b)
                v = NF.rn_old_text(v)
                if v != c.value: changes.append((old, ws.title, c.coordinate, 'M', c.value, v)); c.value = v
    out = os.path.join(OUT, '%s %s.xlsx' % (F['name'], TS))
    wb.save(out)
    # ---------------- checks
    wb2 = openpyxl.load_workbook(out)
    inv = {}
    for ws in wb2.worksheets:
        o = ws.title
        for a, b in ren.items():
            if b.endswith('_') and o.startswith(b): o = a + o[len(b):]
            elif o == b: o = a
        inv[ws.title] = o
    nnum = bad = 0
    for ws in wb2.worksheets:
        old_cells = before[inv[ws.title]]
        for row in ws.iter_rows():
            for c in row:
                if c.value is None or isinstance(c.value, str): continue
                nnum += 1
                if old_cells.get(c.coordinate) != c.value: bad += 1; report.append(('NUMBER CHANGED', old, ws.title, c.coordinate))
        # number cells that disappeared (other than the deleted column)
        for co, v in old_cells.items():
            if isinstance(v, str) or v is None: continue
            if ws[co].value != v and not (old == 3 and ws.title == 'Descriptors' and co.startswith('AJ')):
                bad += 1; report.append(('NUMBER LOST', old, ws.title, co))
    semi, hist, refs = [], [], []
    for ws in wb2.worksheets:
        for row in ws.iter_rows():
            for c in row:
                v = c.value
                if not isinstance(v, str): continue
                if ';' in v and not re.search(r'\[[^\]]*;[^\]]*\]|\$\(', v): semi.append((ws.title, c.coordinate, v[:60]))
                if HIST.search(v): hist.append((ws.title, c.coordinate, v[:80]))
                for m in N.REF_RE.finditer(v):
                    for k, n, l in N.parse(m.group(1), m.group(2)):
                        lim = {'T': 3, 'F': 6, 'TS': 26, 'FS': 8, 'AF': 18}[k]
                        if not 1 <= n <= lim: refs.append((ws.title, c.coordinate, m.group(0)))
    report.append((old, F['new'], os.path.basename(out), 'number cells %d, changed %d' % (nnum, bad),
                   'semicolons left %d' % len(semi), semi[:3], 'history words left %d' % len(hist), hist[:4], 'invalid references %d' % len(refs), refs[:3]))
json.dump(changes, open(os.path.join(OUT, 'Aenderungsliste_AF3_13 %s.json' % TS), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for r in report: print(r)
print('changes', len(changes), '| by kind', {k: sum(1 for c in changes if k in c[3].split()) for k in 'T H S K M'.split()})
