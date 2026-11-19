# -*- coding: utf-8 -*-
"""Check of the descriptor pool of the second data set (37 peptides of Lone et al. 2019).

Recomputes the 1,830 two-dimensional descriptors (Mordred 2D plus RDKit) from the SMILES in
revision/lone/rows.json, applies the cleaning steps stated in Additional file 5 (columns with more than
50 % missing values removed, remaining gaps imputed with the column median, constant columns removed)
and compares the result with the archived matrix revision/lone/desc_pool.csv: column names, column
order and values. The CSV is read with float_precision='round_trip', so the stored decimal values are
compared exactly.

Usage:  python mordred_peptides_check.py <archive_root> <out.json> [<per_descriptor.csv>]
Needs mordredcommunity and rdkit. Verified with mordredcommunity 2.0.7 and 2.0.6, RDKit 2026.09.1 and
2026.03.6, numpy 2.5.3, pandas 3.0.5, Python 3.13.18. Writes nothing into the archive.
"""
import json, os, sys, warnings
from importlib.metadata import version
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from mordred import Calculator, descriptors
from mordred.error import Missing, Error
warnings.filterwarnings('ignore'); RDLogger.DisableLog('rdApp.*')

root, out = sys.argv[1], sys.argv[2]
rows = json.load(open(os.path.join(root, 'revision', 'lone', 'rows.json')))
pool = pd.read_csv(os.path.join(root, 'revision', 'lone', 'desc_pool.csv'), float_precision='round_trip')
assert list(pool['No'].astype(str)) == [r['No'] for r in rows], 'row order of rows.json and desc_pool.csv differs'
pool = pool.drop(columns='No')
mols = [Chem.MolFromSmiles(r['SMILES']) for r in rows]
assert all(m is not None for m in mols)

calc = Calculator(descriptors, ignore_3D=True)
md = calc.pandas(mols, nproc=1, quiet=True)
md.columns = ['Mordred::' + str(d) for d in calc.descriptors]
rd = pd.DataFrame([[f(m) for _, f in Descriptors._descList] for m in mols],
                  columns=['RDKit::' + n for n, _ in Descriptors._descList])
raw = pd.concat([md.reset_index(drop=True), rd], axis=1)
num = raw.apply(lambda col: [np.nan if isinstance(v, (Missing, Error)) else float(v) for v in col]).astype(float)
num = num.replace([np.inf, -np.inf], np.nan)
frac = num.isna().mean()
kept = num.loc[:, frac <= 0.5]
imputed_cols = [c for c in kept.columns if kept[c].isna().any()]
n_imputed = int(kept.isna().sum().sum())
filled = kept.fillna(kept.median())
clean = filled.loc[:, filled.nunique(dropna=False) > 1]

res = dict(versions=dict(mordredcommunity=version('mordredcommunity'), rdkit=Chem.rdBase.rdkitVersion,
                         numpy=np.__version__, pandas=pd.__version__),
           n_peptides=len(mols), raw_descriptors=int(raw.shape[1]),
           removed_more_than_half_missing=int((frac > 0.5).sum()), removed_all_missing=int((frac == 1).sum()),
           removed_constant=int(filled.shape[1] - clean.shape[1]), cleaned=int(clean.shape[1]),
           archived=int(pool.shape[1]), imputed_cells=n_imputed, imputed_columns=imputed_cols,
           imputed_peptides=sorted({rows[i]['No'] for c in imputed_cols for i in np.where(kept[c].isna())[0]}),
           same_columns_same_order=list(clean.columns) == list(pool.columns))
if res['same_columns_same_order']:
    A, B = clean.to_numpy(float), pool.to_numpy(float)
    ad = np.abs(A - B); rel = ad / np.maximum(np.maximum(np.abs(A), np.abs(B)), 1e-300); rel[ad == 0] = 0
    res.update(cells=int(A.size), bit_identical_cells=int((ad == 0).sum()), max_relative_difference=float(rel.max()),
               cells_above_1e_12=int((rel > 1e-12).sum()),
               descriptors_above_1e_9=sorted({clean.columns[j] for i, j in np.argwhere(rel > 1e-9)}))
    if len(sys.argv) > 3:
        pd.DataFrame(dict(descriptor=clean.columns, max_abs_diff=ad.max(axis=0), max_rel_diff=rel.max(axis=0),
                          imputed=[c in imputed_cols for c in clean.columns])).to_csv(sys.argv[3], index=False)
else:
    res.update(only_recomputed=sorted(set(clean.columns) - set(pool.columns)),
               only_archived=sorted(set(pool.columns) - set(clean.columns)))
json.dump(res, open(out, 'w'), indent=1)
for k, v in res.items():
    print('%-34s %s' % (k, v))
