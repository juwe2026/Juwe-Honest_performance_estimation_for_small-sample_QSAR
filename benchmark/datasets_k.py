"""Data sets of the main analysis following the supplementary specification (without control and reference substances).
Assignment of the cells from kontrollen/auswahl_neu.json (mechanical repetition of the selection), exclusions from
kontrollen/entscheidungen.json. Pipeline unchanged from datasets.py and teilmengen.py.
Unchanged data sets (no compound removed, same assay) are taken over from datasets/.
Output: datasets_k/DSxx_{X.csv, y.csv, info.json}
"""
import json, os, shutil, time, sys
import numpy as np, pandas as pd
import importlib.util
spec = importlib.util.spec_from_file_location('ds', 'datasets.py'); ds = importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)

OUT = 'datasets_k'; os.makedirs(OUT, exist_ok=True)
SEL = json.load(open('kontrollen/auswahl_neu.json'))['belegt']
ENT = json.load(open('kontrollen/entscheidungen.json'))
CANDS = {x['assay']: x for x in json.load(open('cand_ranked.json'))}
ALT = {json.load(open(f'datasets/DS{i:02d}_info.json'))['info']['assay']: f'DS{i:02d}' for i in range(1, 25) if os.path.exists(f'datasets/DS{i:02d}_info.json')}


def baue(name, zelle, assay, art, dokument, auswahl=None, seed=None):
    if os.path.exists(f'{OUT}/{name}_info.json'):
        print(name, 'already there'); return
    t0 = time.time()
    aus = set(ENT.get(dokument, {}).get('ausschluss', []))
    rows, unit = ds.klassen(assay, art)
    entfernt = [r[1] for r in rows if r[1] in aus]
    rows = [r for r in rows if r[1] not in aus]
    eltern_n = len(rows)
    if auswahl:
        n, n_minor = auswahl
        y_all = np.array([r[4] for r in rows])
        minor_label = 1 if y_all.sum() <= (y_all == 0).sum() else 0
        i_min = np.where(y_all == minor_label)[0]; i_maj = np.where(y_all != minor_label)[0]
        rng = np.random.default_rng(seed)
        pick = np.sort(np.concatenate([rng.choice(i_min, size=n_minor, replace=False), rng.choice(i_maj, size=n - n_minor, replace=False)]))
        rows = [rows[i] for i in pick]
    y = np.array([r[4] for r in rows]); mols = [r[3] for r in rows]
    X, roh = ds.reinigen(ds.deskriptoren(mols))
    gruppen, cols = ds.gruppieren(X)
    x = CANDS[assay]
    info = dict(ds=name, zelle=list(zelle), dokument=dokument, doi=x.get('doi'), assay=assay,
                ziel=(x['organismus'] or x['ziel']) if art == 'MIC' else (x['ziel'] or x['organismus']),
                einheit=unit, n=len(rows), aktiv=int(y.sum()), nicht_aktiv=int((y == 0).sum()),
                anteil_kleinere=round(min(int(y.sum()), int((y == 0).sum())) / len(y), 4),
                p_roh=roh, p_sauber=X.shape[1], gruppen=len(gruppen), redundanz_G_p=round(len(gruppen) / X.shape[1], 4),
                entfernt=entfernt, eltern_n=eltern_n if auswahl else None, seed=seed, sekunden=round(time.time() - t0, 1))
    X.to_csv(f'{OUT}/{name}_X.csv', index=False)
    pd.DataFrame(dict(molregno=[r[0] for r in rows], chembl_id=[r[1] for r in rows], smiles=[r[2] for r in rows], y=y,
                      relation=[r[5] for r in rows], wert=[r[6] for r in rows])).to_csv(f'{OUT}/{name}_y.csv', index=False)
    json.dump(dict(info=info, gruppen=[[cols[j] for j in g] for g in gruppen]), open(f'{OUT}/{name}_info.json', 'w'), indent=1)
    print(name, info['zelle'], info['dokument'], 'n', info['n'], f"{info['aktiv']}/{info['nicht_aktiv']}", 'p', info['p_sauber'],
          'G', info['gruppen'], 'G/p', info['redundanz_G_p'], 'removed', len(entfernt), f"{info['sekunden']}s", flush=True)


if __name__ == '__main__':
    zellen = sorted(SEL.items(), key=lambda kv: eval(kv[0]))
    for i, (k, v) in enumerate(zellen, 1):
        name = f'DS{i:02d}'; zelle = eval(k)
        alt = ALT.get(v['assay'])
        if alt == name and not ENT.get(v['dokument'], {}).get('ausschluss'):
            for suf in ('_X.csv', '_y.csv', '_info.json'):
                shutil.copy(f'datasets/{name}{suf}', f'{OUT}/{name}{suf}')
            print(name, 'taken over unchanged'); continue
        baue(name, zelle, v['assay'], zelle[0], v['dokument'])
    ELTERN = {'MIC': dict(dokument='CHEMBL3407479', assay='CHEMBL3414284'), 'Enzym': dict(dokument='CHEMBL5360622', assay='CHEMBL5363956')}
    for name, art, n, n_minor, dsnum in [('DS21', 'MIC', 25, 8, 21), ('DS22', 'MIC', 60, 18, 22), ('DS23', 'Enzym', 25, 8, 23), ('DS24', 'Enzym', 60, 18, 24)]:
        e = ELTERN[art]
        baue(name, ('Teilmenge', art, f'n={n}'), e['assay'], art, e['dokument'], auswahl=(n, n_minor), seed=[93, dsnum, 9 + 200])
