"""Builds the four subsets of the size variation (DS21 to DS24).
From each parent data set two stratified draws, n = 25 and n = 60, with 30 % in the minor class.
No models, no permutations.
"""
import json, os, time
import numpy as np, pandas as pd
import importlib.util

spec = importlib.util.spec_from_file_location('ds', 'datasets.py')
ds = importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)

ELTERN = json.load(open('elternpools.json'))
PLAN = [('DS21', 'MIC', 25, 8, 21), ('DS22', 'MIC', 60, 18, 22),
        ('DS23', 'Enzym', 25, 8, 23), ('DS24', 'Enzym', 60, 18, 24)]
OUT = 'datasets'

cache = {}
for name, art, n, n_minor, dsnum in PLAN:
    if os.path.exists(f'{OUT}/{name}_info.json'):
        print(name, 'already there'); continue
    t0 = time.time()
    e = ELTERN[art]
    if art not in cache:
        rows, unit = ds.klassen(e['assay'], art)
        cache[art] = (rows, unit)
        print(f'parent set {art}: {len(rows)} compounds read', flush=True)
    rows, unit = cache[art]
    y_all = np.array([r[4] for r in rows])
    # the minor class of the parent set determines which class is drawn to n_minor
    minor_label = 1 if y_all.sum() <= (y_all == 0).sum() else 0
    idx_minor = np.where(y_all == minor_label)[0]
    idx_major = np.where(y_all != minor_label)[0]
    rng = np.random.default_rng([93, dsnum, 9])
    pick_minor = rng.choice(idx_minor, size=n_minor, replace=False)
    pick_major = rng.choice(idx_major, size=n - n_minor, replace=False)
    pick = np.sort(np.concatenate([pick_minor, pick_major]))
    sub = [rows[i] for i in pick]
    y = np.array([r[4] for r in sub]); mols = [r[3] for r in sub]
    df = ds.deskriptoren(mols)
    X, roh = ds.reinigen(df)
    gruppen, cols = ds.gruppieren(X)
    info = dict(ds=name, zelle=['Teilmenge', art, f'n={n}'], dokument=e['dokument'], doi=e['doi'],
                assay=e['assay'], ziel=e['ziel'], einheit=unit, n=len(sub),
                aktiv=int(y.sum()), nicht_aktiv=int((y == 0).sum()),
                anteil_kleinere=round(min(int(y.sum()), int((y == 0).sum())) / len(y), 4),
                p_roh=roh, p_sauber=X.shape[1], gruppen=len(gruppen),
                redundanz_G_p=round(len(gruppen) / X.shape[1], 4),
                eltern_n=len(rows), seed=[93, dsnum, 9], sekunden=round(time.time() - t0, 1))
    X.to_csv(f'{OUT}/{name}_X.csv', index=False)
    pd.DataFrame(dict(molregno=[r[0] for r in sub], chembl_id=[r[1] for r in sub],
                      smiles=[r[2] for r in sub], y=y, relation=[r[5] for r in sub],
                      wert=[r[6] for r in sub])).to_csv(f'{OUT}/{name}_y.csv', index=False)
    json.dump(dict(info=info, gruppen=[[cols[j] for j in g] for g in gruppen]),
              open(f'{OUT}/{name}_info.json', 'w'), indent=1)
    print(name, art, 'n', info['n'], f"{info['aktiv']}/{info['nicht_aktiv']}", 'p', info['p_sauber'],
          'G', info['gruppen'], 'G/p', info['redundanz_G_p'], f"{info['sekunden']}s", flush=True)
