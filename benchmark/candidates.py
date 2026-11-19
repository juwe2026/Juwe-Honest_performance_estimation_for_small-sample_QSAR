"""Candidate list of the ChEMBL benchmark (Revision 2) following the pre-specification, draft 3.
No models, no filters, no permutations. Only selection key figures.
Usage: python3 candidates.py chembl37_subset.db out_prefix
Memory-saving: activities are streamed per assay.
"""
import sys, sqlite3, hashlib, math, statistics, json
from collections import defaultdict, Counter
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')

db = sqlite3.connect(sys.argv[1]); out = sys.argv[2]
q = lambda s, a=(): db.execute(s, a).fetchall()

# ---------- protein classes ----------
pc = {r[0]: (r[1], r[2]) for r in q("select protein_class_id,parent_id,pref_name from protein_classification")}
def lineage(cid):
    names = []
    while cid is not None and cid in pc:
        p, n = pc[cid]; names.append(n); cid = p
    return names
tid_classes = defaultdict(set)
for tid, comp in q("select tc.tid, cc.protein_class_id from target_components tc join component_class cc on cc.component_id=tc.component_id"):
    tid_classes[int(tid)].add(tuple(lineage(comp)))
def prot_kind(tid):
    kinds = set()
    for lin in tid_classes.get(tid, ()):
        s = [x.lower() for x in lin]
        if 'enzyme' in s: kinds.add('Enzym')
        if 'membrane receptor' in s or 'ion channel' in s or 'nuclear receptor' in s: kinds.add('Rezeptor')
    return kinds

targets = {int(r[0]): r[1:] for r in q("select tid,target_type,pref_name,organism,chembl_id from target_dictionary")}
NONBACT = ('plasmodium', 'trypanosoma', 'leishmania', 'toxoplasma', 'schistosoma', 'virus', 'homo sapiens', 'mus musculus',
           'rattus', 'giardia', 'entamoeba', 'trichomonas', 'brugia', 'haemonchus', 'caenorhabditis', 'babesia', 'cryptosporidium',
           'onchocerca', 'eimeria', 'naegleria', 'acanthamoeba', 'fasciola', 'echinococcus', 'artemia', 'aedes', 'anopheles', 'culex')
docs = {int(r[0]): r for r in q("select doc_id,journal,year,volume,first_page,pubmed_id,doi,chembl_id,title,doc_type,src_id,authors from docs where doc_type='PUBLICATION'")}

# ---------- admissible assays ----------
db.execute("drop table if exists elig")
db.execute("create table elig(assay_id integer primary key, ttype text, pk text)")
buf = []
for aid, did, tid, conf in db.execute("select assay_id,doc_id,tid,confidence_score from assays"):
    if did is None or int(did) not in docs or not tid: continue
    tid = int(tid)
    if tid not in targets: continue
    ttype, tname, torg, tchembl = targets[tid]
    if ttype == 'ORGANISM':
        org = (torg or tname or '').lower()
        if not org or any(x in org for x in NONBACT): continue
        buf.append((int(aid), 'ORG', ''))
    elif ttype == 'CELL-LINE':
        buf.append((int(aid), 'CELL', ''))
    elif ttype == 'SINGLE PROTEIN' and conf not in (None, '') and int(conf) >= 8:
        k = prot_kind(tid)
        if k == {'Enzym'}: buf.append((int(aid), 'SP', 'Enzym'))
        elif k == {'Rezeptor'}: buf.append((int(aid), 'SP', 'Rezeptor'))
db.executemany("insert into elig values (?,?,?)", buf); db.commit(); del buf
db.execute("create index if not exists i_act_assay on activities(assay_id)"); db.commit()
print('admissible assays', q("select ttype,pk,count(*) from elig group by 1,2"), flush=True)

def kind_of(stype, ttype, pk):
    if ttype == 'ORG' and stype == 'MIC': return 'MIC'
    if ttype == 'CELL' and stype in ('IC50', 'GI50', 'CC50'): return 'Zytotox'
    if ttype == 'SP' and stype == 'IC50' and pk == 'Enzym': return 'Enzym'
    if ttype == 'SP' and stype == 'Ki' and pk == 'Rezeptor': return 'Rezeptor'
    return None

# ---------- structures ----------
def parent_of(m):
    r = q("select parent_molregno from molecule_hierarchy where molregno=?", (str(m),))
    return int(r[0][0]) if r and r[0][0] else m
METALS = set('Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Cs Ba La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po'.split())
_mc = {}
def compound_ok(m):
    if m in _mc:
        return Chem.MolFromSmiles(_mc[m]) if _mc[m] else None
    mol = None
    md = q("select structure_type,inorganic_flag,polymer_flag from molecule_dictionary where molregno=?", (str(m),))
    if md and md[0][0] == 'MOL' and md[0][1] != '1' and md[0][2] != '1':
        r = q("select canonical_smiles from compound_structures where molregno=?", (str(m),))
        smi = r[0][0] if r else None
        if smi and '.' not in smi:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None and (mol.GetNumHeavyAtoms() < 6 or any(at.GetSymbol() in METALS for at in mol.GetAtoms())):
                mol = None
    _mc[m] = Chem.MolToSmiles(mol) if mol is not None else None
    return mol

FALLBACK = {'MIC': ('ugml', 8.0), 'Enzym': ('p', 6.0), 'Rezeptor': ('p', 7.0), 'Zytotox': ('p', 5.0)}
mfp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def iter_groups():
    cur = db.execute("""select a.assay_id,a.molregno,a.standard_relation,a.standard_value,a.standard_units,a.standard_type,e.ttype,e.pk
        from activities a join elig e on e.assay_id=a.assay_id
        where a.src_id='1' and a.data_validity_comment is null and a.potential_duplicate='0'
        and a.standard_value is not null and a.molregno is not null""")
    # avoid grouping in SQLite: collecting per assay only works when sorted, hence an intermediate table
    return cur

db.execute("drop table if exists act_f")
db.execute("""create table act_f as select cast(a.assay_id as integer) aid, cast(a.molregno as integer) molregno, a.standard_relation rel,
        cast(a.standard_value as real) val, a.standard_units unit, a.standard_type stype, e.ttype ttype, e.pk pk
        from activities a join elig e on e.assay_id=a.assay_id
        where a.src_id='1' and a.data_validity_comment is null and a.potential_duplicate='0'
        and a.standard_value is not null and a.molregno is not null""")
db.execute("create index i_actf on act_f(aid, stype)"); db.commit()
grp = q("select aid, stype, count(distinct molregno) from act_f group by aid, stype having count(distinct molregno) >= 30")
print('groups with >=30 compounds', len(grp), flush=True)

cands = []
for gi, (aid, stype, nraw) in enumerate(grp):
    if gi % 5000 == 0: print(gi, len(cands), flush=True)
    rs = q("select molregno,rel,val,unit,ttype,pk from act_f where aid=? and stype=?", (aid, stype))
    kind = kind_of(stype, rs[0][4], rs[0][5])
    if kind is None: continue
    units = {r[3] for r in rs}
    if len(units) != 1: continue
    unit = units.pop()
    if kind == 'MIC' and unit not in ('ug.mL-1', 'nM'): continue
    if kind != 'MIC' and unit != 'nM': continue
    per = defaultdict(list)
    for m, rel, v, _, _, _ in rs:
        per[parent_of(m)].append((rel, v))
    if not (30 <= len(per) <= 1000): continue
    mols = {}; vals = {}
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
    n = len(vals)
    if n < 30: continue
    mode, thr = FALLBACK[kind]
    act = nonact = undecided = 0
    for m, (rel, v) in vals.items():
        if v is None or v <= 0:
            undecided += 1; continue
        if mode == 'ugml':
            x = v if unit == 'ug.mL-1' else v * 1e-9 * Descriptors.MolWt(mols[m]) * 1e3
            ok_act = x <= thr
        else:
            ok_act = (9 - math.log10(v)) >= thr
        if rel == '=':
            act += ok_act; nonact += (not ok_act)
        elif rel == '>':
            if not ok_act: nonact += 1
            else: undecided += 1
        else:
            if ok_act: act += 1
            else: undecided += 1
    scaf = Counter(); fps = []
    for m, mol in mols.items():
        try:
            g = MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(mol))
            scaf[Chem.MolToSmiles(g)] += 1
        except Exception:
            scaf['?'] += 1
        fps.append(mfp.GetFingerprint(mol))
    top = scaf.most_common(1)[0][1] / len(mols)
    tot = 0.0; cnt = 0
    for i in range(len(fps) - 1):
        s = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]); tot += sum(s); cnt += len(s)
    mts = tot / cnt
    div = 'kongener' if top >= 0.5 else ('divers' if (top < 0.25 and mts < 0.30) else 'dazwischen')
    a = q("select assay_id,doc_id,description,assay_type,assay_organism,tid,confidence_score,chembl_id from assays where assay_id=?", (str(aid),))[0]
    d = docs[int(a[1])]; t = targets[int(a[5])]
    key = hashlib.sha256(f"QSAR-R2-{d[7]}{a[7]}".encode()).hexdigest()
    ndec = act + nonact; minor = min(act, nonact)
    cands.append(dict(hash=key, art=kind, standard_type=stype, einheit=unit, assay=a[7], dokument=d[7], journal=d[1], jahr=d[2],
                      pmid=d[5], doi=d[6], titel=d[8], ziel=t[1], ziel_chembl=t[3], organismus=t[2], n_verbindungen=n,
                      n_unentschieden_vorl=undecided, aktiv_vorl=act, nichtaktiv_vorl=nonact,
                      anteil_kleinere_klasse_vorl=round(minor / ndec, 3) if ndec else None,
                      gerueste=len(scaf), groesstes_geruest_anteil=round(top, 3), mittl_tanimoto=round(mts, 3), diversitaet=div,
                      beschreibung=(a[2] or '')[:300]))
print('candidates after structure check', len(cands), Counter(c['art'] for c in cands), flush=True)
json.dump(cands, open(out + '_all.json', 'w'), ensure_ascii=False)
