"""Mechanical repetition of the selection with the supplementary rule (control and reference substances).

Procedure exactly as select.py (hash order, per publication the first matching assay in a free band, one data set per
publication, per target, second pass with repeated target). New: before an assay is checked, the
foreign substances are removed and size, classes and diversity are recomputed. No model is computed.

Foreign substances:
- for publications already read and decided, the decision from kontrollen/entscheidungen.json
- otherwise a mechanical search criterion: compound with a ChEMBL preferred name and highest Tanimoto similarity
  to the unnamed compounds of the assay below 0.60. If an assay cleaned in this way is selected or its status changes as a result,
  the article is read and the decision is entered in entscheidungen.json.
"""
import json, sqlite3, math, statistics, sys
from collections import defaultdict, Counter
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')

db = sqlite3.connect('chembl37_subset.db')
q = lambda s, a=(): db.execute(s, a).fetchall()
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
FALLBACK = {'MIC': ('ugml', 8.0), 'Enzym': ('p', 6.0), 'Rezeptor': ('p', 7.0), 'Zytotox': ('p', 5.0)}
METALS = set('Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Cs Ba La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po'.split())
ENT = json.load(open('kontrollen/entscheidungen.json'))   # document -> {"gelesen": true, "ausschluss": [chembl_id, ...]}   (gelesen = read, ausschluss = excluded)
SUCH_GRENZE = 0.60

PAR = {int(m): int(p) for m, p in db.execute("select molregno,parent_molregno from molecule_hierarchy") if p}
NAMED = {int(m) for (m,) in db.execute("select molregno from molecule_dictionary where pref_name is not null and pref_name!=''")}
def parent_of(m):
    return PAR.get(int(m), int(m))
N_VOLL = [0]

_mc = {}
def compound_ok(m):
    if m in _mc: return _mc[m]
    mol = None
    md = q("select structure_type,inorganic_flag,polymer_flag from molecule_dictionary where molregno=?", (str(m),))
    if md and md[0][0] == 'MOL' and md[0][1] != '1' and md[0][2] != '1':
        r = q("select canonical_smiles from compound_structures where molregno=?", (str(m),))
        smi = r[0][0] if r else None
        if smi and '.' not in smi:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None and (mol.GetNumHeavyAtoms() < 6 or any(at.GetSymbol() in METALS for at in mol.GetAtoms())):
                mol = None
    _mc[m] = mol
    return mol

_info = {}
def molinfo(m):
    if m not in _info:
        r = q("select chembl_id, pref_name from molecule_dictionary where molregno=?", (str(m),))
        _info[m] = (r[0][0], r[0][1] or '') if r else ('', '')
    return _info[m]

_aid = {}
def assay_num(chembl):
    if chembl not in _aid:
        _aid[chembl] = int(q("select assay_id from assays where chembl_id=?", (chembl,))[0][0])
    return _aid[chembl]

_cache = {}
def bewerte(x):
    """Assay x (entry from cand_ranked.json) after removal of the foreign substances."""
    key = x['assay']
    if key in _cache: return _cache[key]
    rs = q("select molregno,rel,val,unit from act_f where aid=? and stype=?", (assay_num(x['assay']), x['standard_type']))
    if x['dokument'] not in ENT and not any(parent_of(m) in NAMED for m, _, _, _ in rs):
        # no named compound: key figures unchanged as in the original candidate list
        a_, n_ = x['aktiv_vorl'], x['nichtaktiv_vorl']
        r = dict(n=x['n_verbindungen'], n_vorher=x['n_verbindungen'], aktiv=a_, nicht=n_, unentschieden=x['n_unentschieden_vorl'],
                 anteil=x['anteil_kleinere_klasse_vorl'], minor=min(a_, n_), top=x['groesstes_geruest_anteil'],
                 mts=x['mittl_tanimoto'], div=x['diversitaet'], quelle='unveraendert', entfernt=[])
        _cache[key] = r
        return r
    # An assay that still has more than 60 compounds after removal of all named compounds cannot fill a cell.
    # Then the original assessment suffices (saves the costly similarity computation for large assays).
    k_named = len({parent_of(m) for m, _, _, _ in rs if parent_of(m) in NAMED})
    if x['dokument'] not in ENT and x['n_verbindungen'] - k_named > 60:
        a_, n_ = x['aktiv_vorl'], x['nichtaktiv_vorl']
        r = dict(n=x['n_verbindungen'], n_vorher=x['n_verbindungen'], aktiv=a_, nicht=n_, unentschieden=x['n_unentschieden_vorl'],
                 anteil=x['anteil_kleinere_klasse_vorl'], minor=min(a_, n_), top=x['groesstes_geruest_anteil'],
                 mts=x['mittl_tanimoto'], div=x['diversitaet'], quelle='zu gross auch ohne benannte Verbindungen', entfernt=[])
        _cache[key] = r
        return r
    N_VOLL[0] += 1
    if N_VOLL[0] % 50 == 0: print('fully assessed', N_VOLL[0], flush=True)
    unit = x['einheit']; kind = x['art']
    per = defaultdict(list)
    for m, rel, v, u in rs:
        per[parent_of(m)].append((rel, v))
    mols, vals = {}, {}
    for m, lst in per.items():
        mol = compound_ok(m)
        if mol is None: continue
        eq = [v for rel, v in lst if rel == '=' and v and v > 0]
        if eq:
            vals[m] = ('=', 10 ** statistics.median([math.log10(v) for v in eq]))
        else:
            rels = {rel for rel, _ in lst}
            if rels <= {'>', '>='}: vals[m] = ('>', max(v for _, v in lst))
            elif rels <= {'<', '<='}: vals[m] = ('<', min(v for _, v in lst))
            else: continue
        mols[m] = mol
    ids = list(mols)
    fps = {m: GEN.GetFingerprint(mols[m]) for m in ids}
    # determine the foreign substances
    ent = ENT.get(x['dokument'])
    weg = []
    if ent and ent.get('gelesen'):
        aus = set(ent.get('ausschluss', []))
        weg = [m for m in ids if molinfo(m)[0] in aus]
        quelle = 'gelesen'
    else:
        # comparison only with the unnamed compounds, so that several similar reference substances
        # (e.g. several azoles) do not mark each other as a series
        unbenannt = [fps[k] for k in ids if not molinfo(k)[1]]
        for m in ids:
            if molinfo(m)[1]:
                if not unbenannt or max(DataStructs.BulkTanimotoSimilarity(fps[m], unbenannt)) < SUCH_GRENZE:
                    weg.append(m)
        quelle = 'Suchkriterium'
    rest = [m for m in ids if m not in weg]
    mode, thr = FALLBACK[kind]
    act = nonact = und = 0
    for m in rest:
        rel, v = vals[m]
        if v is None or v <= 0: und += 1; continue
        if mode == 'ugml':
            xv = v if unit == 'ug.mL-1' else v * 1e-9 * Descriptors.MolWt(mols[m]) * 1e3
            ok = xv <= thr
        else:
            ok = (9 - math.log10(v)) >= thr
        if rel == '=': act += ok; nonact += (not ok)
        elif rel == '>':
            if not ok: nonact += 1
            else: und += 1
        else:
            if ok: act += 1
            else: und += 1
    scaf = Counter()
    for m in rest:
        try: scaf[Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(mols[m])))] += 1
        except Exception: scaf['?'] += 1
    top = scaf.most_common(1)[0][1] / len(rest) if rest else 0
    tot = cnt = 0
    fl = [fps[m] for m in rest]
    for i in range(len(fl) - 1):
        s = DataStructs.BulkTanimotoSimilarity(fl[i], fl[i + 1:]); tot += sum(s); cnt += len(s)
    mts = tot / cnt if cnt else 0
    div = 'kongener' if top >= 0.5 else ('divers' if (top < 0.25 and mts < 0.30) else 'dazwischen')
    n = len(rest)
    nd = act + nonact
    r = dict(n=n, n_vorher=len(ids), aktiv=act, nicht=nonact, unentschieden=und,
             anteil=round(min(act, nonact) / nd, 4) if nd else None, minor=min(act, nonact),
             top=round(top, 4), mts=round(mts, 5), div=div, quelle=quelle,
             entfernt=[(molinfo(m)[0], molinfo(m)[1]) for m in weg])
    _cache[key] = r
    return r


ARTS = ['MIC', 'Enzym', 'Rezeptor', 'Zytotox']; DIVS = ['kongener', 'divers']; BANDS = ['A 25-35 %', 'B >35-50 %']
def band_of(f):
    if f is None or f < 0.25 or f > 0.50: return None
    return BANDS[0] if f <= 0.35 else BANDS[1]
def ziel(x):
    return (x['organismus'] or x['ziel']) if x['art'] == 'MIC' else (x['ziel'] or x['organismus'])

def pruefe(x):
    """(cell or None, reason) after cleaning."""
    r = bewerte(x)
    if not (30 <= r['n'] <= 60): return None, f"Groesse {r['n']} ausserhalb 30-60"
    if r['div'] == 'dazwischen': return None, f"Diversitaet dazwischen (top {r['top']}, Tanimoto {r['mts']})"
    if r['minor'] < 10: return None, f"Klassengrenze verfehlt ({r['anteil']}, kleinere Klasse {r['minor']})"
    b = band_of(r['anteil'])
    if b is None: return None, f"Klassengrenze verfehlt ({r['anteil']})"
    return (x['art'], r['div'], b), ''

def lauf(cands, verdikte, bis_alle=True):
    bydoc = defaultdict(list)
    for x in cands:
        if x['grund'].startswith('Groesse >= 100'): continue
        bydoc[x['dokument']].append(x)
    for d in bydoc: bydoc[d].sort(key=lambda y: y['hash'])
    docs = sorted(bydoc, key=lambda d: bydoc[d][0]['hash'])
    belegt, used_doc, used_t, prot, lesen = {}, set(), set(), [], []
    def offen(): return [(a, dv, b) for a in ARTS for dv in DIVS for b in BANDS if (a, dv, b) not in belegt]
    for d in docs:
        if not offen(): break
        v = verdikte.get(d, {}).get('art')
        if v in ('review', 'methoden', 'unlesbar', 'keine_serie', 'vorselektiert'):
            continue
        # check only publications that have any assay of a cell type that is still open
        if not any(x['art'] in {z[0] for z in offen()} for x in bydoc[d]): continue
        gew = None
        for i, x in enumerate(bydoc[d], 1):
            if x['art'] not in {z[0] for z in offen()}: continue
            zelle, grund = pruefe(x)
            if zelle is None: prot.append((d, x['assay'], grund)); continue
            if ziel(x) in used_t: prot.append((d, x['assay'], f'Ziel schon belegt: {ziel(x)}')); continue
            if zelle in belegt: prot.append((d, x['assay'], f'Band schon besetzt: {zelle}')); continue
            gew = (x, zelle, i); break
        if gew:
            x, zelle, i = gew
            if v is None:
                lesen.append((d, x['assay'], zelle)); prot.append((d, x['assay'], f'MUSS GELESEN WERDEN fuer {zelle}'))
                # the cell stays reserved until the article has been read: mark it as provisionally occupied
                belegt[zelle] = dict(x, _vorlaeufig=True, _rang=i); used_doc.add(d); used_t.add(ziel(x)); continue
            belegt[zelle] = dict(x, _rang=i); used_doc.add(d); used_t.add(ziel(x))
            prot.append((d, x['assay'], f'GEWAEHLT {zelle}'))
    # second pass
    tc = Counter(ziel(x) for x in belegt.values())
    for d in docs:
        if not offen(): break
        v = verdikte.get(d, {}).get('art')
        if v in ('review', 'methoden', 'unlesbar', 'keine_serie', 'vorselektiert') or d in used_doc: continue
        for i, x in enumerate(bydoc[d], 1):
            if x['art'] not in {z[0] for z in offen()}: continue
            zelle, grund = pruefe(x)
            if zelle is None or zelle in belegt or tc[ziel(x)] >= 2: continue
            if v is None:
                lesen.append((d, x['assay'], zelle)); belegt[zelle] = dict(x, _vorlaeufig=True, _rang=i, _durchgang=2)
            else:
                belegt[zelle] = dict(x, _rang=i, _durchgang=2); prot.append((d, x['assay'], f'GEWAEHLT im 2. Durchgang {zelle}'))
            used_doc.add(d); tc[ziel(x)] += 1; break
    return belegt, prot, lesen, offen()


if __name__ == '__main__':
    cands = json.load(open('cand_ranked.json'))
    verdikte = json.load(open('verdikte.json'))
    belegt, prot, lesen, offen = lauf(cands, verdikte)
    for z in sorted(belegt):
        x = belegt[z]; r = bewerte(x)
        print(z, x['dokument'], x['assay'], ziel(x)[:30], 'n', r['n'], f"{r['aktiv']}/{r['not']}", r['anteil'], r['div'],
              'removed', [e[1] or e[0] for e in r['entfernt']], 'PROVISIONAL' if x.get('_vorlaeufig') else '')
    print('open:', offen)
    print('to be read:', lesen)
    json.dump(dict(belegt={str(k): dict(dokument=v['dokument'], assay=v['assay'], vorlaeufig=v.get('_vorlaeufig', False),
                                         bewertung=bewerte(v)) for k, v in belegt.items()},
                   protokoll=prot, lesen=lesen, offen=[list(o) for o in offen]),
              open('kontrollen/auswahl_neu.json', 'w'), indent=1, ensure_ascii=False)
