# -*- coding: utf-8 -*-
"""Approved numbering in citation order (EW 06.11.2026 13:34 Berlin, 'Nummerierung_Zitierreihenfolge_VORSCHLAG
06.11.2026 1325', N2 with N1). Maps old (working) number -> new number. renumber_ms.py recomputes the order from the
text and stops if it differs from these maps."""

T = {1: 1, 5: 2, 6: 3}
F = {1: 1, 2: 2, 3: 3, 6: 4, 7: 5, 8: 6}
TS_ORDER = [1, 2, 3, 4, 23, 5, 6, 7, 8, 9, 10, 19, 11, 13, 15, 24, 25, 26, 17, 12, 20, 21, 14, 18, 22, 16]
FS_ORDER = [2, 3, 4, 7, 8, 1, 5, 6]
AF_ORDER = [1, 2, 6, 3, 5, 4, 13, 15, 14, 7, 11, 16, 8, 9, 10, 12, 17, 18]
TS = {o: i + 1 for i, o in enumerate(TS_ORDER)}
FS = {o: i + 1 for i, o in enumerate(FS_ORDER)}
AF = {o: i + 1 for i, o in enumerate(AF_ORDER)}
MAPS = {'T': T, 'F': F, 'TS': TS, 'FS': FS, 'AF': AF}
# items cited only in Additional file 1, placed after the items cited in the main text
ONLY_AF1 = {'TS': [16], 'FS': []}
# old main tables that moved to Additional file 1 (working numbers), for references in Additional files 1 and 2
MOVED_MAIN_TABLES = {2: 24, 3: 25, 4: 26}
MOVED_MAIN_FIGURES = {4: 7, 5: 8}

def rn_old_text(t):
    """Text in working numbers -> new numbers, including the former main Tables 2 to 4 (now supplementary)."""
    import re, nummern as N
    T_ext = dict(T); T_ext.update({k: 9000 + k for k in MOVED_MAIN_TABLES})
    maps = dict(MAPS); maps['T'] = T_ext
    t = N.sub_text(t, maps)
    return re.sub(r'\bTable 900(\d)\b', lambda m: 'Table S%d' % TS[MOVED_MAIN_TABLES[int(m.group(1))]], t)
