# -*- coding: utf-8 -*-
"""Renumbering of tables, figures and additional files in citation order (decision EW 06.11.2026 06:55 and 13:16).

Kinds: 'T' main tables, 'F' main figures, 'TS' supplementary tables, 'FS' supplementary figures, 'AF' additional files.
A reference is 'Table 1', 'Tables S1 and S2', 'Fig. 8d', 'Figs. S5 and S6', 'Table S1 to S22', 'Additional file 7',
'Additional files 2, 3, 6 and 8 to 12'. Lists are mapped item by item. A list without panel letters is written
sorted, with runs of three or more as 'a to b' and the last item joined by 'and'. Every list whose form changes
beyond the substitution of numbers is reported."""
import re

NUM = r'S?\d+[a-z]?(?![.,]\d)'          # '7, 0.41' is not a list, '1,286' is not a number list
SEP = r'(?:, | and | to )'
REF_RE = re.compile(r'\b(Tables?|Figs?\.|Additional files?) (' + NUM + r'(?:' + SEP + NUM + r')*)\b')

def kind_of(word, item):
    s = item.startswith('S')
    if word.startswith('Table'): return 'TS' if s else 'T'
    if word.startswith('Fig'): return 'FS' if s else 'F'
    return 'AF'

def parse(word, lst):
    """-> list of (kind, number, panel letter), ranges expanded. Raises on mixed kinds."""
    toks = re.split(r'(, | and | to )', lst)
    items, rng = [], False
    for t in toks:
        if t == ' to ': rng = True; continue
        if t in (', ', ' and '): continue
        m = re.fullmatch(r'(S?)(\d+)([a-z]?)', t); assert m, (word, lst)
        k = kind_of(word, t)
        if rng:
            pk, pn, pl = items[-1]; assert pk == k and not pl and not m.group(3), (word, lst)
            items += [(k, n, '') for n in range(pn + 1, int(m.group(2)) + 1)]
            rng = False
        else:
            items.append((k, int(m.group(2)), m.group(3)))
    assert len({k for k, n, l in items}) == 1, ('mixed list', word, lst)
    return items

def fmt_list(kind, nums):
    """Sorted numbers -> 'S1, S3 and S5 to S8' (prefix S for supplementary kinds)."""
    pre = 'S' if kind in ('TS', 'FS') else ''
    runs, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        if j - i >= 2: runs.append('%s%d to %s%d' % (pre, nums[i], pre, nums[j]))
        else: runs += ['%s%d' % (pre, n) for n in nums[i:j + 1]]
        i = j + 1
    return runs[0] if len(runs) == 1 else ', '.join(runs[:-1]) + ' and ' + runs[-1]

def rewrite(word, lst, maps):
    """-> (new word, new list text, changed form?)."""
    items = parse(word, lst)
    kind = items[0][0]
    mp = maps.get(kind)
    if mp is None: return word, lst, False
    new = [(k, mp[n], l) for k, n, l in items]
    pre = 'S' if kind in ('TS', 'FS') else ''
    plain = re.sub(r'(S?)(\d+)', lambda m: m.group(1) + str(mp[int(m.group(2))]), lst)
    if len(items) == 1 or any(l for k, n, l in items):
        return word, plain, False
    nums = sorted(n for k, n, l in new)
    assert len(set(nums)) == len(nums)
    out = fmt_list(kind, nums)
    return word, out, out != plain

def sub_text(t, maps, log=None, where=''):
    def f(m):
        w, l, changed = rewrite(m.group(1), m.group(2), maps)
        if changed and log is not None: log.append((where, m.group(0), w + ' ' + l))
        return w + ' ' + l
    return REF_RE.sub(f, t)

def first_order(texts, kinds=('T', 'F', 'TS', 'FS', 'AF'), skip_own_label=True):
    """texts: paragraph texts in document order. -> {kind: [old numbers in order of first citation]}."""
    order = {k: [] for k in kinds}
    for t in texts:
        own = re.match(r'(Fig\. \d+|Table \d+)\b', t) if skip_own_label else None
        for m in REF_RE.finditer(t):
            if own and m.start() == 0: continue
            for k, n, l in parse(m.group(1), m.group(2)):
                if k in order and n not in order[k]: order[k].append(n)
    return order
