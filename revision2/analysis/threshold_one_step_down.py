# -*- coding: utf-8 -*-
"""Active compounds of the monoterpenoids when the activity threshold is moved one dilution step down
(main text, Methods, item G1 of 16.11.2026).

The activity classes of the main analysis count a compound as active when the mean MIC of its liquid-phase
replicates is below 2048 µg/mL. One step down means a mean below 1024 µg/mL. For comparison the script also gives
the classes under the median of the replicates at 512 µg/mL or below. It reads the replicate values from the sheet
MIC_replicates_31 of Additional file 2.

Usage:  python threshold_one_step_down.py <Additional file 2 (.xlsx)>
Nothing is written. Running time a few seconds."""
import sys, collections, statistics as st
import openpyxl

wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
rows = list(wb['MIC_replicates_31'].iter_rows(values_only=True))
head, rows = rows[0], rows[1:]
print('columns:', head)
vals = collections.defaultdict(list)
for compound, strain, phase, replicate, mic in rows:
    if str(phase).startswith('liquid'):
        vals[(strain, compound)].append(float(mic))
rules = [('mean below 2048 (main analysis)', lambda v: st.mean(v) < 2048),
         ('mean below 1024 (one step down)', lambda v: st.mean(v) < 1024),
         ('median at or below 512', lambda v: st.median(v) <= 512)]
for strain in sorted({s for s, _ in vals}):
    comps = sorted(c for s, c in vals if s == strain)
    print('%s (%d compounds)' % (strain, len(comps)))
    for name, f in rules:
        act = [c for c in comps if f(vals[(strain, c)])]
        print('   %-32s %2d  %s' % (name, len(act), ', '.join(act) if len(act) < 8 else ''))
    for c in comps:
        m = st.mean(vals[(strain, c)])
        if 512 < m < 1024:
            print('   between 512 and 1024: %s, mean %.2f, median %g, %d values' % (c, m, st.median(vals[(strain, c)]), len(vals[(strain, c)])))
