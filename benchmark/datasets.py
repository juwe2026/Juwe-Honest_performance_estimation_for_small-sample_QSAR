"""Builds the 16 data sets of the ChEMBL benchmark following the frozen specification.
No models, no permutations. Only compounds, classes, descriptors, grouping.
Output: datasets/DSxx_{X.csv, y.csv, info.json}
"""
import json, sqlite3, math, statistics, os, sys, time
from collections import defaultdict, Counter
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from scipy.stats import spearmanr
from mordred import Calculator, descriptors as mdesc
RDLogger.DisableLog('rdApp.*')

DB = 'chembl37_subset.db'
OUT = 'datasets'
os.makedirs(OUT, exist_ok=True)
db = sqlite3.connect(DB)
q = lambda s, a=(): db.execute(s, a).fetchall()

FALLBACK = {'MIC': ('ugml', 8.0), 'Enzym': ('p', 6.0), 'Rezeptor': ('p', 7.0), 'Zytotox': ('p', 5.0)}
METALS = set('Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Cs Ba La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po'.split())

def parent_of(m):
    r = q("select parent_molregno from molecule_hierarchy where molregno=?", (str(m),))
    return int(r[0][0]) if r and r[0][0] else m

def mol_of(m):
    md = q("select structure_type,inorganic_flag,polymer_flag,chembl_id from molecule_dictionary where molregno=?", (str(m),))
    if not md or md[0][0] != 'MOL' or md[0][1] == '1' or md[0][2] == '1': return None, None
    r = q("select canonical_smiles from compound_structures where molregno=?", (str(m),))
    smi = r[0][0] if r else None
    if not smi or '.' in smi: return None, None
    mol = Chem.MolFromSmiles(smi)
    if mol is None or mol.GetNumHeavyAtoms() < 6: return None, None
    if any(a.GetSymbol() in METALS for a in mol.GetAtoms()): return None, None
    return mol, md[0][3]

def klassen(assay_chembl, art):
    aid = q("select assay_id from assays where chembl_id=?", (assay_chembl,))[0][0]
    rows = q("""select molregno, standard_relation, standard_value, standard_units from activities
                where assay_id=? and src_id='1' and data_validity_comment is null and potential_duplicate='0'
                and standard_value is not null and molregno is not null""", (str(aid),))
    units = {r[3] for r in rows}
    assert len(units) == 1, f'inconsistent unit: {units}'
    unit = units.pop()
    per = defaultdict(list)
    for m, rel, val, _ in rows:
        per[parent_of(int(m))].append((rel, float(val)))
    mode, thr = FALLBACK[art]
    out = []
    for m, lst in sorted(per.items()):
        mol, cid = mol_of(m)
        if mol is None: continue
        eq = [v for rel, v in lst if rel == '=' and v > 0]
        if eq:
            rel, v = '=', 10 ** statistics.median([math.log10(x) for x in eq])
        else:
            rels = {r for r, _ in lst}
            if rels <= {'>', '>='}: rel, v = '>', max(v for _, v in lst)
            elif rels <= {'<', '<='}: rel, v = '<', min(v for _, v in lst)
            else: continue
        if v is None or v <= 0: continue
        if mode == 'ugml':
            x = v if unit == 'ug.mL-1' else v * 1e-9 * Descriptors.MolWt(mol) * 1e3
            akt = x <= thr
        else:
            akt = (9 - math.log10(v)) >= thr
        if rel == '=':
            y = 1 if akt else 0
        elif rel == '>':
            if akt: continue          # undecidable
            y = 0
        else:
            if not akt: continue
            y = 1
        out.append((m, cid, Chem.MolToSmiles(mol), mol, y, rel, v, unit))
    return out, unit

CALC = Calculator(mdesc, ignore_3D=True)
RD = [d for d in Descriptors._descList]

def deskriptoren(mols):
    md = CALC.pandas(mols, quiet=True, nproc=1)
    md = md.astype(object)
    import pandas as pd
    rd = pd.DataFrame([{n: f(m) for n, f in RD} for m in mols])
    rd.columns = ['RD_' + c for c in rd.columns]
    return pd.concat([md.reset_index(drop=True), rd.reset_index(drop=True)], axis=1)

def reinigen(df):
    import pandas as pd
    roh = df.shape[1]
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.loc[:, df.isna().mean() <= 0.5]
    df = df.fillna(df.median(numeric_only=True))
    df = df.loc[:, df.std(ddof=0) > 1e-12]
    df = df.loc[:, ~df.columns.duplicated()]
    return df, roh

def gruppieren(X):
    cols = list(X.columns); A = X.values.astype(float)
    P = np.corrcoef(A, rowvar=False)
    S = spearmanr(A).statistic
    if np.ndim(S) == 0: S = np.array([[1.0]])
    R = (np.abs(np.nan_to_num(P)) > 0.7) | (np.abs(np.nan_to_num(S)) > 0.7)
    np.fill_diagonal(R, False)
    var = A.var(axis=0, ddof=0)
    offen = set(range(len(cols))); gruppen = []
    while offen:
        grade = {i: sum(1 for j in offen if j != i and R[i, j]) for i in offen}
        mx = max(grade.values())
        kand = [i for i in offen if grade[i] == mx]
        hub = max(kand, key=lambda i: (var[i], -i))
        mit = {j for j in offen if j != hub and R[hub, j]}
        gruppen.append([hub] + sorted(mit))
        offen -= ({hub} | mit)
    return gruppen, cols

if __name__ == '__main__':
    import pandas as pd
    import importlib.util
    spec = importlib.util.spec_from_file_location('sel', 'select.py'); sel = importlib.util.module_from_spec(spec); spec.loader.exec_module(sel)
    c = json.load(open('cand_ranked.json')); v = json.load(open('verdikte.json'))
    belegt, _, _, _ = sel.run(c, v)
    nur = sys.argv[1:] or None
    zus = []
    for i, (k, x) in enumerate(sorted(belegt.items()), 1):
        ds = f'DS{i:02d}'
        if nur and ds not in nur: continue
        if os.path.exists(f'{OUT}/{ds}_info.json'): print(ds, 'already there'); continue
        t0 = time.time()
        rows, unit = klassen(x['assay'], x['art'])
        y = np.array([r[4] for r in rows]); mols = [r[3] for r in rows]
        df = deskriptoren(mols)
        X, roh = reinigen(df)
        gruppen, cols = gruppieren(X)
        info = dict(ds=ds, zelle=list(k), dokument=x['dokument'], doi=x.get('doi'), assay=x['assay'],
                    ziel=(x['organismus'] or x['ziel']) if x['art'] == 'MIC' else (x['ziel'] or x['organismus']),
                    einheit=unit, n=len(rows), aktiv=int(y.sum()), nicht_aktiv=int((y == 0).sum()),
                    anteil_kleinere=round(min(int(y.sum()), int((y == 0).sum())) / len(y), 4),
                    p_roh=roh, p_sauber=X.shape[1], gruppen=len(gruppen),
                    redundanz_G_p=round(len(gruppen) / X.shape[1], 4), sekunden=round(time.time() - t0, 1))
        X.to_csv(f'{OUT}/{ds}_X.csv', index=False)
        pd.DataFrame(dict(molregno=[r[0] for r in rows], chembl_id=[r[1] for r in rows],
                          smiles=[r[2] for r in rows], y=y, relation=[r[5] for r in rows],
                          wert=[r[6] for r in rows])).to_csv(f'{OUT}/{ds}_y.csv', index=False)
        json.dump(dict(info=info, gruppen=[[cols[j] for j in g] for g in gruppen]), open(f'{OUT}/{ds}_info.json', 'w'), indent=1)
        zus.append(info)
        print(ds, info['zelle'], 'n', info['n'], f"{info['aktiv']}/{info['nicht_aktiv']}",
              'p', info['p_sauber'], 'G', info['gruppen'], 'G/p', info['redundanz_G_p'], f"{info['sekunden']}s", flush=True)
    if zus: json.dump(zus, open(f'{OUT}/_zusammenfassung_{int(time.time())}.json', 'w'), indent=1)
