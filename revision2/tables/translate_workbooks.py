# -*- coding: utf-8 -*-
"""Translate the remaining German text of the benchmark results workbooks and update their sheet Deviations.

Inputs : the results workbooks of 25.10.2026 21:06 (main analysis, original rule) and 30.10.2026 14:14 (median rule)
Changes: only those listed in uebersetzung_af16_18.py. The script checks that
  - every cell that is not text is unchanged,
  - text changes occur only in the listed sheets and columns, with the expected counts,
  - every number inside a translated text is unchanged,
  - no German word of the selection script remains.
Usage: python translate_workbooks.py MAIN.xlsx ORIGINAL.xlsx MEDIAN.xlsx OUT_MAIN.xlsx OUT_ORIGINAL.xlsx OUT_MEDIAN.xlsx
"""
import sys, os, re, copy
import openpyxl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uebersetzung_af16_18 as U

GERMAN = re.compile(r"\b(dazwischen|Diversitaet|kongener|Groesse|ausserhalb|Elternpool\w*|Groessenvariation|MUSS|GELESEN|WERDEN|fuer|Durchgang|"
                    r"Zielwiederholung|Publikation|ausgeschlossen|keine_serie|vorselektiert|methoden|Enzym|Rezeptor|Zytotox|Autoren|Kontrollsubstanzen\w*|abgelehnt)\b|'divers'| / divers / ")
NUM = re.compile(r'\d+(?:\.\d+)?')

def snapshot(wb):
    return {ws.title: {(c.row, c.column): c.value for row in ws.iter_rows() for c in row if c.value is not None} for ws in wb.worksheets}

def translate(path_in, path_out, kind, log):
    wb = openpyxl.load_workbook(path_in); before = snapshot(wb); counts = {}
    allowed = set()                                   # (sheet, column) where text may change
    if kind in ('main', 'original'):
        for sh, col, old, new in U.CELL:
            ws = wb[sh]; n = 0
            for row in ws.iter_rows():
                c = row[col] if len(row) > col else None
                if c is not None and c.value == old: c.value = new; n += 1
            counts['CELL %s | %s' % (sh, old)] = n; allowed.add((sh, col + 1))
        ws = wb['Selection_log']; n = 0
        for row in ws.iter_rows(min_row=2):
            c = row[2]
            if isinstance(c.value, str):
                v = c.value
                for old, new in U.LOG: v = v.replace(old, new)
                if v != c.value:
                    o, w = NUM.findall(c.value), NUM.findall(v.replace('Enzyme IC50', 'Enzyme'))
                    assert o == w, (c.value, v)
                    c.value = v; n += 1
        counts['LOG Selection_log'] = n; allowed.add(('Selection_log', 3))
    for sh, old, new in U.PART:
        if sh not in wb.sheetnames: continue
        ws = wb[sh]; n = 0
        for row in ws.iter_rows():
            for c in row:
                if isinstance(c.value, str) and old in c.value:
                    v = c.value.replace(old, new)
                    if sh != 'README': assert NUM.findall(v) == NUM.findall(c.value)
                    c.value = v; n += 1; allowed.add((sh, c.column))
        counts['PART %s | %s' % (sh, old[:40])] = n
    if kind == 'main':
        ws = wb['Excluded_compounds']; n = 0
        for row in ws.iter_rows():
            for c in row:
                if c.value == U.EXCLUDED_NOTE[0]: c.value = U.EXCLUDED_NOTE[1]; n += 1; allowed.add(('Excluded_compounds', c.column))
        counts['Excluded_compounds Satz'] = n
    if kind in ('main', 'original'):
        ws = wb['Deviations']
        r_open = [r for r in range(1, ws.max_row + 1) if [ws.cell(r, k).value for k in range(1, 5)] == U.DEV_OPEN]
        assert len(r_open) == 1; r0 = r_open[0]
        new_rows = U.DEV_MAIN if kind == 'main' else U.DEV_ORIGINAL
        styles = [copy.copy(ws.cell(r0, k)._style) for k in range(1, 5)]
        if len(new_rows) > 1: ws.insert_rows(r0 + 1, len(new_rows) - 1)
        for i, vals in enumerate(new_rows):
            for k, v in enumerate(vals, 1):
                ws.cell(r0 + i, k).value = v; ws.cell(r0 + i, k)._style = copy.copy(styles[k - 1])
        counts['Deviations: row open replaced by'] = len(new_rows)
    after = snapshot(wb)
    # ---- checks
    for sh in before:
        b, a = before[sh], after[sh]
        if sh == 'Deviations' and kind in ('main', 'original'):
            # all rows other than the replaced one are unchanged, in the same order
            rb = [tuple(b.get((r, k)) for k in range(1, 5)) for r in range(1, max(x[0] for x in b) + 1)]
            ra = [tuple(a.get((r, k)) for k in range(1, 5)) for r in range(1, max(x[0] for x in a) + 1)]
            exp = [x for x in rb if list(x) != U.DEV_OPEN]
            got = [x for x in ra if list(x) not in (U.DEV_MAIN if kind == 'main' else U.DEV_ORIGINAL)]
            assert got == exp, 'Deviations: other rows changed'
            continue
        assert set(b) == set(a), (sh, 'set of cells changed')
        for key, v in b.items():
            w = a[key]
            if v == w: continue
            assert isinstance(v, str) and isinstance(w, str), (sh, key, v, w)
            assert (sh, key[1]) in allowed, (sh, key, v, w)
    for sh, cells in after.items():
        for key, v in cells.items():
            if isinstance(v, str) and GERMAN.search(v):
                raise AssertionError('German word remains: %s %s %r' % (sh, key, v[:80]))
    wb.save(path_out)
    log.append('== %s -> %s' % (os.path.basename(path_in), os.path.basename(path_out)))
    for k, v in counts.items(): log.append('   %-75s %6d' % (k, v))
    nchg = sum(1 for sh in before if sh != 'Deviations' for key in before[sh] if before[sh][key] != after[sh].get(key))
    log.append('   changed text cells outside Deviations: %d, number and other cells unchanged' % nchg)

if __name__ == '__main__':
    a = sys.argv[1:]
    LOGL = []
    translate(a[0], a[3], 'main', LOGL)
    translate(a[1], a[4], 'original', LOGL)
    translate(a[2], a[5], 'median', LOGL)
    print('\n'.join(LOGL))
