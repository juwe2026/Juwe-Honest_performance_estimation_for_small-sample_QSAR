"""Data sets of the sensitivity analysis median rule (main analysis, DRAFT of the addendum of 25.10.2026).
Same compounds and descriptor matrices as datasets_k/, only the classes are assigned again following the median rule:
boundary = median of the recorded activity values of the series (per compound summarised as in datasets.py, censored values with their bound).
active = value ('=') strictly smaller than the median. Equal to the median = not active.
Censored values as in the pre-specification: '>' with bound >= median not active, otherwise excluded. '<' with bound < median active, otherwise excluded.
If a compound is excluded, the data set is rebuilt without it with the pipeline of datasets.py (cleaning, grouping).
Output: datasets_m/DSxx_{X.csv, y.csv, info.json}. y.csv additionally contains the column y_ersatzregel (class under the fallback rule).
"""
import json, os, shutil, importlib.util
import numpy as np, pandas as pd
spec = importlib.util.spec_from_file_location('ds', 'datasets.py'); ds = importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)
OUT = 'datasets_m'; os.makedirs(OUT, exist_ok=True)
proto = {}
for i in range(1, 17):
    d = f'DS{i:02d}'
    y = pd.read_csv(f'datasets_k/{d}_y.csv'); info = json.load(open(f'datasets_k/{d}_info.json'))
    v = y.wert.values.astype(float); rel = y.relation.values; med = float(np.median(v))
    lv = np.log10(v); lm = np.log10(med); TOL = 1e-9          # comparison on the logarithmic scale with a tolerance against rounding noise
    unter = lv < lm - TOL; ueber_gleich = lv >= lm - TOL
    ym = np.full(len(y), -1)
    ym[(rel == '=') & unter] = 1; ym[(rel == '=') & ueber_gleich] = 0
    ym[(rel == '>') & ueber_gleich] = 0; ym[(rel == '<') & unter] = 1
    weg = list(y.chembl_id[ym < 0])
    y2 = y.copy(); y2['y_ersatzregel'] = y2['y']; y2['y'] = ym
    if not weg:
        shutil.copy(f'datasets_k/{d}_X.csv', f'{OUT}/{d}_X.csv')
        y2.to_csv(f'{OUT}/{d}_y.csv', index=False)
        I = info['info']
    else:
        keep = ym >= 0; y2 = y2[keep].reset_index(drop=True)
        mols = [ds.Chem.MolFromSmiles(s) for s in y2.smiles]
        X, roh = ds.reinigen(ds.deskriptoren(mols)); gruppen, cols = ds.gruppieren(X)
        X.to_csv(f'{OUT}/{d}_X.csv', index=False); y2.to_csv(f'{OUT}/{d}_y.csv', index=False)
        I = dict(info['info']); I.update(n=len(y2), p_roh=roh, p_sauber=X.shape[1], gruppen=len(gruppen), redundanz_G_p=round(len(gruppen) / X.shape[1], 4))
        info = dict(info, gruppen=[[cols[j] for j in g] for g in gruppen])
    I = dict(I, aktiv=int((y2.y == 1).sum()), nicht_aktiv=int((y2.y == 0).sum()), anteil_kleinere=round(min((y2.y == 1).mean(), (y2.y == 0).mean()), 4),
             median_grenze=med, median_ausgeschlossen=weg, klassen_geaendert=int((y2.y != y2.y_ersatzregel).sum()), regel='Median-Regel')
    json.dump(dict(info, info=I), open(f'{OUT}/{d}_info.json', 'w'), indent=1)
    proto[d] = dict(n=I['n'], aktiv=I['aktiv'], median=med, ausgeschlossen=weg, geaendert=I['klassen_geaendert'], p=I['p_sauber'], G=I['gruppen'])
    print(d, proto[d], flush=True)
json.dump(proto, open(f'{OUT}/protokoll.json', 'w'), indent=1)
