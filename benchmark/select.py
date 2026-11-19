"""Selection following the pre-specification, draft 5.
Publications in hash order. Within a publication, assays in hash order.
The first assay that meets the class boundary AND falls into a free band is taken.
One data set per publication, one data set per organism/target/cell line.
Bands: A = 25 to 35 %, B = above 35 to 50 %.
Verdicts from reading are in verdikte.json: document -> {"art": "original|review|methoden|unlesbar", "schwelle": ...}
(art = type of article: original, review, methods paper, unreadable. schwelle = threshold of the authors.)
"""
import json, sys
from collections import defaultdict

ARTS = ['MIC', 'Enzym', 'Rezeptor', 'Zytotox']
DIVS = ['kongener', 'divers']
BANDS = ['A 25-35 %', 'B >35-50 %']

def band_of(frac):
    if frac is None: return None
    if frac < 0.25 or frac > 0.50: return None
    return BANDS[0] if frac <= 0.35 else BANDS[1]

def qualifies(x):
    a, n = x['aktiv_vorl'], x['nichtaktiv_vorl']
    tot = a + n
    if not tot: return None, None
    mi = min(a, n)
    if mi < 10: return None, None
    frac = mi / tot
    return band_of(frac), frac

def run(cands, verdikte):
    # group publications, rank = rank of the first assay
    bydoc = defaultdict(list)
    for x in cands:
        if x['grund']: continue          # size or diversity not admissible
        bydoc[x['dokument']].append(x)
    for d in bydoc: bydoc[d].sort(key=lambda y: y['hash'])
    docs = sorted(bydoc, key=lambda d: bydoc[d][0]['hash'])

    belegt = {}          # (type, diversity, band) -> data set
    used_doc, used_target = set(), set()
    protokoll = []
    for d in docs:
        v = verdikte.get(d, {}).get('art')
        if v in ('review', 'methoden', 'unlesbar', 'keine_serie', 'vorselektiert'):
            protokoll.append((d, None, f'Publikation ausgeschlossen: {v}')); continue
        if v is None:
            continue      # not read yet, cannot be assigned yet
        if d in used_doc: continue
        gewaehlt = None
        for i, x in enumerate(bydoc[d], 1):
            b, frac = qualifies(x)
            tgt = (x['organismus'] or x['ziel']) if x['art'] == 'MIC' else (x['ziel'] or x['organismus'])
            key = (x['art'], x['diversitaet'], b)
            if b is None:
                protokoll.append((d, x['assay'], f'Klassengrenze verfehlt ({x["anteil_kleinere_klasse_vorl"]})')); continue
            if tgt in used_target:
                protokoll.append((d, x['assay'], f'Ziel schon belegt: {tgt}')); continue
            if key in belegt:
                protokoll.append((d, x['assay'], f'Band schon besetzt: {key}')); continue
            gewaehlt = (x, b, i); break
        if gewaehlt:
            x, b, i = gewaehlt
            belegt[(x['art'], x['diversitaet'], b)] = dict(x, _band=b, _assay_rang=i)
            used_doc.add(d); used_target.add((x['organismus'] or x['ziel']) if x['art'] == 'MIC' else (x['ziel'] or x['organismus']))
            protokoll.append((d, x['assay'], f'GEWAEHLT {x["art"]} / {x["diversitaet"]} / {b}'))
    # second pass: the target may repeat if a cell would otherwise stay empty (at most twice per target)
    from collections import Counter
    tgt_count = Counter()
    for x in belegt.values():
        tgt_count[(x['organismus'] or x['ziel']) if x['art'] == 'MIC' else (x['ziel'] or x['organismus'])] += 1
    for d in docs:
        if not offene_zellen(belegt): break
        v = verdikte.get(d, {}).get('art')
        if v != 'original' or d in used_doc: continue
        for i, x in enumerate(bydoc[d], 1):
            b, frac = qualifies(x)
            if b is None: continue
            key = (x['art'], x['diversitaet'], b)
            if key in belegt: continue
            tgt = (x['organismus'] or x['ziel']) if x['art'] == 'MIC' else (x['ziel'] or x['organismus'])
            if tgt_count[tgt] >= 2: continue
            belegt[key] = dict(x, _band=b, _assay_rang=i, _durchgang=2)
            used_doc.add(d); tgt_count[tgt] += 1
            protokoll.append((d, x['assay'], f'GEWAEHLT im 2. Durchgang (Zielwiederholung) {key}'))
            break
    return belegt, protokoll, bydoc, docs

def offene_zellen(belegt):
    out = []
    for a in ARTS:
        for dv in DIVS:
            for b in BANDS:
                if (a, dv, b) not in belegt: out.append((a, dv, b))
    return out

if __name__ == '__main__':
    cands = json.load(open('cand_ranked.json'))
    try: verdikte = json.load(open('verdikte.json'))
    except FileNotFoundError: verdikte = {}
    belegt, prot, bydoc, docs = run(cands, verdikte)
    print('occupied:', len(belegt))
    for k, v in belegt.items():
        print('  ', k, v['assay'], ((v['organismus'] or v['ziel']) if v['art']=='MIC' else (v['ziel'] or v['organismus']))[:34], v['n_verbindungen'],
              f"{v['aktiv_vorl']}/{v['nichtaktiv_vorl']}", v['anteil_kleinere_klasse_vorl'], 'assay rank', v['_assay_rang'])
    print('open:', len(offene_zellen(belegt)))
    for z in offene_zellen(belegt): print('  ', z)
