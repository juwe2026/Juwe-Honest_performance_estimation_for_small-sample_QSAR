# -*- coding: utf-8 -*-
"""Check of the main figures against the BMC figure requirements (submission guidelines of the Journal of
Cheminformatics, 'Figure size and resolution', as quoted by EW on 06.11.2026): width 85 or 170 mm, height of figure
and legend at most 225 mm, about 300 dpi at final size, lines wider than 0.25 pt, text legible.

For each vector PDF: native size, all text spans with their font size, all stroked paths with their line width
(PyMuPDF). Line widths are those set in the drawing, scaled by the current transformation (PyMuPDF reports them in
page units). Results are given at native size and at the scale factors that follow from the limits."""
import sys, json, collections
import pymupdf
MM = 25.4 / 72
out = {}
for name in sys.argv[1:]:
    doc = pymupdf.open(name); pg = doc[0]
    w, h = pg.rect.width * MM, pg.rect.height * MM
    sizes = collections.Counter()
    texts = collections.defaultdict(list)
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if s['text'].strip():
                    sz = round(s['size'], 2); sizes[sz] += len(s['text'].strip()); texts[sz].append(s['text'].strip())
    lw = collections.Counter()
    for d in pg.get_drawings():
        if d.get('type') in ('s', 'fs') and d.get('width') is not None:
            lw[round(d['width'], 3)] += 1
    out[name] = dict(w_mm=round(w, 1), h_mm=round(h, 1), font_sizes={str(k): v for k, v in sorted(sizes.items())},
                     smallest_font=min(sizes), smallest_font_examples=texts[min(sizes)][:6],
                     line_widths={str(k): v for k, v in sorted(lw.items())}, thinnest_line=min(lw) if lw else None)
    print('%s  %.1f x %.1f mm | font sizes (pt: characters) %s | smallest %s %s | line widths (pt: paths) %s' % (
        name, w, h, dict(sorted(sizes.items())), min(sizes), texts[min(sizes)][:4], dict(sorted(lw.items()))))
json.dump(out, open('abb_pruefung.json', 'w'), indent=1, ensure_ascii=False)
