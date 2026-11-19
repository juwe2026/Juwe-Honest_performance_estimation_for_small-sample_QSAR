# Revision (G2-14, EW 03.10.2026): ionisation and tautomerism, Route B.
# (1) cis-pinonic acid and geranic acid as carboxylate: TPSA change (RDKit, Ertl) and Route B Q2 (1 and 3 LV) with otherwise
#     unchanged descriptors. (2) beta-thujaplicin: the two tropolone tautomers, TPSA, aromatic rings, Phenol_flag SMARTS.
# Deterministic. Needs QSAR_ROUTEA_XLSX and QSAR_ROUTEB_XLSX.
# Archive-relative paths. _ROOT is the archive root, found from the location of this file. _REV is revision/,
# _OUT holds the outputs of the revision scripts (default results/revision, or QSAR_REV_OUT), and _WORK holds
# intermediate files (default results/recalc, or QSAR_RECALC_DIR).
import os as _osp
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_REV = _osp.path.join(_ROOT, 'revision')
_OUT = _osp.environ.get('QSAR_REV_OUT', _osp.path.join(_ROOT, 'results', 'revision'))
_WORK = _osp.environ.get('QSAR_RECALC_DIR', _osp.path.join(_ROOT, 'results', 'recalc'))
_osp.makedirs(_WORK, exist_ok=True)
exec(open((_REV + '/common.py')).read())
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
json.dump(out, open((_OUT + '/ionisation_check.json'), 'w'), indent=1, default=float)
print({n: round(v, 2) for n, v in dT.items()}, 'max |dQ2| 1LV %.4f' % max(abs(v['diff']) for k, v in out['routeB_q2'].items() if k.endswith('_1LV')))
print(out['thujaplicin'])
