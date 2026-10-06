# Revision (G2-14, EW 03.10.2026): Ionisation und Tautomerie, Route B.
# (1) cis-Pinonsaeure und Geraniumsaeure als Carboxylat: TPSA-Aenderung (RDKit, Ertl) und Route-B-Q2 (1 und 3 LV) mit sonst
#     unveraenderten Deskriptoren. (2) beta-Thujaplicin: die beiden Tropolon-Tautomere, TPSA, aromatische Ringe, Phenol_flag-SMARTS.
# Deterministisch. Braucht QSAR_ROUTEA_XLSX und QSAR_ROUTEB_XLSX.
exec(open('/home/claude/rev/common.py').read())
import core, json, openpyxl
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rd
wb = openpyxl.load_workbook(os.environ['QSAR_ROUTEB_XLSX'], read_only=True, data_only=True)
rows = list(wb['Descriptors'].iter_rows(values_only=True)); smi = {r[0]: r[1] for r in rows[1:] if r[0]}
out = {'acids': {}, 'routeB_q2': {}, 'thujaplicin': []}
for n in ('cis-Pinonic acid', 'Geranic acid'):
    m = Chem.MolFromSmiles(smi[n]); mi = Chem.MolFromSmiles(Chem.MolToSmiles(m).replace('C(=O)O', 'C(=O)[O-]'))
    out['acids'][n] = dict(smiles=smi[n], tpsa=rd.CalcTPSA(m), tpsa_ion=rd.CalcTPSA(mi))
dT = {n: v['tpsa_ion'] - v['tpsa'] for n, v in out['acids'].items()}
X12 = core.X_of(core.RED12).copy(); j = core.RED12.index('TPSA'); Xi = X12.copy()
for i, n in enumerate(core.COMPOUNDS):
    if n in dT: Xi[i, j] += dT[n]
S = {'HI': 'H.I.', 'SA': 'S.A.', 'SPneu': 'S.Pneu', 'SPyo': 'S.Pyo', 'PA': 'P.A.'}
for L in (1, 3):
    for k in KEYS:
        y = core.y_of(S[k]); a = q2(y, loo_pls(X12, y, L)); b = q2(y, loo_pls(Xi, y, L))
        out['routeB_q2']['%s_%dLV' % (k, L)] = dict(neutral=a, carboxylate=b, diff=b - a)
phen = Chem.MolFromSmarts('[c][OX2H1]')
for s in ('CC(C)c1cccc(=O)c(O)c1', 'CC(C)c1cccc(O)c(=O)c1'):
    t = Chem.MolFromSmiles(s)
    out['thujaplicin'].append(dict(smiles=s, tpsa=rd.CalcTPSA(t), aromatic_rings=rd.CalcNumAromaticRings(t), phenol_flag=int(t.HasSubstructMatch(phen))))
json.dump(out, open('/home/claude/rev/out/ionisation_check.json', 'w'), indent=1, default=float)
print({n: round(v, 2) for n, v in dT.items()}, 'max |dQ2| 1LV %.4f' % max(abs(v['diff']) for k, v in out['routeB_q2'].items() if k.endswith('_1LV')))
print(out['thujaplicin'])
