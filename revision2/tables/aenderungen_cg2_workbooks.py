# -*- coding: utf-8 -*-
"""Corrections to the published Additional files 4, 12, 14, 15 and 17 after the review report of 15 November 2026
(items G12 and G13, approved by EW on 16 November 2026).

G12  P. aeruginosa, which has three active compounds, was called an "artefact" or an "overfitting artefact" in
     Additional files 4, 14 and 15. The main text and Additional file 2 call it "descriptive only". The six text
     cells and eight table cells below now say so. No number changes.
G13  In the sheet Unstable_predictions of Additional files 12 and 17 the column Value(s) holds the ChEMBL relation
     and value as text, for example "=7000". Because the text begins with "=", openpyxl stored 106 and 216 of these
     cells as formulas without a cached value, so that spreadsheet programs compute them and programs that read the
     file show them as empty. They are now stored as text with the same content.

The files are edited in their XML, cell by cell. Every other part of each file stays byte for byte as it was.
The script checks with openpyxl that every cell other than the listed ones is unchanged and that every changed
formula cell now holds exactly the formula text it held before, as a text cell.

Usage: python aenderungen_cg2_workbooks.py IN_DIR OUT_DIR TIMESTAMP
IN_DIR holds the published Additional files 4, 12, 14, 15 and 17 (names beginning with Additional_file_<n>_).
OUT_DIR receives them with the time stamp TIMESTAMP in place of the old one."""
import sys, os, re, glob, zipfile
from xml.sax.saxutils import escape
import openpyxl

TEXT = {   # file number -> list of (sheet, cell, old text, new text)
    4: [('Conclusions', 'A15', '5) P. aeruginosa (3 actives) is an overfitting artefact (AUC = 1.00) and is reported descriptively only.',
                               '5) P. aeruginosa has only three active compounds and is reported descriptively only (AUC = 1.00).'),
        ('Comparison_3_routes', 'G6', 'Artefact (3 actives): descriptive only, no p-value interpreted',
                                      'Three actives: descriptive only, no p-value interpreted')],
    14: [('Info', 'A19', 'P. aeruginosa has three active compounds, is an artefact and is reported descriptively only.',
                         'P. aeruginosa has three active compounds and is reported descriptively only.')] +
        [('Top15_RouteA1_EffSize_vs_VIP', c, 'artefact', 'descriptive') for c in ('E6', 'G6', 'K6', 'M6')],
    15: [('Info', 'A16', 'P. aeruginosa (3 actives) is an artefact and is reported descriptively only.',
                         'P. aeruginosa (3 actives) is reported descriptively only.')] +
        [('Top15_RouteA2_EffSize_vs_VIP', c, 'artefact', 'descriptive') for c in ('D6', 'F6', 'I6', 'K6')],
}
FORMULA = {12: ('Unstable_predictions', 'E', 106), 17: ('Unstable_predictions', 'E', 216)}

SRC, OUT, TS = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)


def sheet_part(z, name):
    wbx = z.read('xl/workbook.xml').decode(); rels = z.read('xl/_rels/workbook.xml.rels').decode()
    rid = re.search(r'<sheet [^>]*name="%s"[^>]*r:id="(rId\d+)"' % re.escape(name), wbx).group(1)
    m = re.search(r'<Relationship [^>]*Target="([^"]+)"[^>]*Id="%s"' % rid, rels) or \
        re.search(r'<Relationship [^>]*Id="%s"[^>]*Target="([^"]+)"' % rid, rels)
    t = m.group(1)
    return t.lstrip('/') if t.startswith('/') else 'xl/' + t


def snapshot(path):
    wb = openpyxl.load_workbook(path)
    return {ws.title: {c.coordinate: (c.value, c.data_type) for row in ws.iter_rows() for c in row if c.value is not None}
            for ws in wb.worksheets}


for num in (4, 12, 14, 15, 17):
    src = glob.glob(os.path.join(SRC, 'Additional_file_%d_*.xlsx' % num))
    assert len(src) == 1, (num, src); src = src[0]
    z = zipfile.ZipFile(src); parts = {}
    expected = {}                              # (sheet, cell) -> (value, data_type) after the change
    for sh, cell, old, new in TEXT.get(num, []):
        p = sheet_part(z, sh); x = parts.get(p) or z.read(p).decode('utf-8')
        pat = re.compile(r'(<c r="%s"[^>]*t="inlineStr"><is><t>)%s(</t></is></c>)' % (cell, re.escape(escape(old))))
        x, n = pat.subn(lambda m: m.group(1) + escape(new) + m.group(2), x); assert n == 1, (num, sh, cell)
        parts[p] = x; expected[(sh, cell)] = (new, 's')
    if num in FORMULA:
        sh, col, count = FORMULA[num]; p = sheet_part(z, sh); x = z.read(p).decode('utf-8')
        assert x.count('<f>') == count, (num, x.count('<f>'))
        pat = re.compile(r'<c r="(%s\d+)"( s="\d+")?><f>([^<]*)</f><v></v></c>' % col)
        found = pat.findall(x); assert len(found) == count and all(re.fullmatch(r'\d+(\.\d+)?', f) for _, _, f in found)
        x = pat.sub(lambda m: '<c r="%s"%s t="inlineStr"><is><t>=%s</t></is></c>' % (m.group(1), m.group(2) or '', m.group(3)), x)
        assert '<f>' not in x
        parts[p] = x
        for cell, _, f in found: expected[(sh, cell)] = ('=' + f, 's')
    name = re.sub(r'\d\d\.\d\d\.\d{4} \d{4}\.xlsx$', TS + '.xlsx', os.path.basename(src)); assert name != os.path.basename(src)
    out = os.path.join(OUT, name)
    with zipfile.ZipFile(out, 'w') as zo:
        for it in z.infolist():
            data = parts[it.filename].encode('utf-8') if it.filename in parts else z.read(it.filename)
            zo.writestr(it, data)
    # ---------------- checks
    a, b = snapshot(src), snapshot(out)
    assert list(a) == list(b)
    nchg = 0
    for sh in a:
        keys = set(a[sh]) | set(b[sh])
        for k in keys:
            if (sh, k) in expected:
                assert b[sh][k] == expected[(sh, k)], (num, sh, k, b[sh][k])
                if num in FORMULA and a[sh][k][1] == 'f': assert a[sh][k][0] == b[sh][k][0]
                nchg += 1
            else:
                assert a[sh].get(k) == b[sh].get(k), (num, sh, k)
    zs, zo = zipfile.ZipFile(src), zipfile.ZipFile(out)
    same = [n for n in zs.namelist() if n not in parts and zs.read(n) == zo.read(n)]
    print('Additional file %d: %s -> %s, cells changed %d of %d expected, other parts byte-identical %d of %d'
          % (num, os.path.basename(src), name, nchg, len(expected), len(same), len(zs.namelist()) - len(parts)))
