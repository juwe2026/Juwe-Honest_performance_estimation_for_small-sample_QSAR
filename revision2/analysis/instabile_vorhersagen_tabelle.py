"""Table DS03, DS11, DS15 (main analysis): all compounds of the assay, control substances and outliers marked.
The cell texts of the working workbook are in German (internal working document, English content in Additional file 12)."""
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
import json, sys, math, statistics, importlib.util
from collections import Counter, defaultdict
import numpy as np, pandas as pd
sys.path.insert(0, (_ROOT + '/benchmark'))
import os; os.chdir((_ROOT + '/benchmark'))
if not os.path.exists('chembl37_subset.db'):
    sys.exit('benchmark/chembl37_subset.db is needed (built from the ChEMBL 37 release, see benchmark/README_benchmark.md)')
spec = importlib.util.spec_from_file_location('ds', 'datasets.py'); ds = importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

OUT = sys.argv[1]
V = json.load(open('ergebnisse_k/vorhersagen.json'))
CR = json.load(open('kontrollen/compound_records.json'))
MOD = {'app_A1': 'A1 scheinbar', 'app_A1u': 'A1u scheinbar', 'app_A2': 'A2 scheinbar', 'A1': 'A1 genestet', 'A1u': 'A1u genestet',
       'A2': 'A2 genestet', 'ohneAuswahl': 'ohne Auswahl', 'ohneFilter': 'ohne Filter'}
GRUND = {
 'CHEMBL2018886': 'Einzige Verbindung mit diesem Grundgerüst in der Serie (1 von 45), mit 36 Schweratomen größer als der Rest (Median 29). Ein Deskriptor ist bei den übrigen Verbindungen fast konstant, bei dieser Verbindung aber deutlich anders. Die Autoskalierung vergrößert diesen Unterschied auf |z| = {z}. Dadurch weit außerhalb liegende Vorhersagen in fünf Modellen.',
 'CHEMBL248514': 'Gehört zum Hauptgerüst der Serie (31 von 34). Ein Deskriptor ist bei den übrigen Verbindungen fast konstant, bei dieser Verbindung aber deutlich anders (|z| = {z}). Betroffen ist nur das Kontrollmodell ohne Filter, alle Modelle mit Auswahl sind unauffällig.',
 'CHEMBL603616': 'Gehört zum Hauptgerüst der Serie (33 von 39). Die Hexadecylkette ist kein Unterscheidungsmerkmal, 14 Verbindungen der Serie tragen sie. Ein Deskriptor ist bei den übrigen Verbindungen fast konstant, bei dieser Verbindung aber deutlich anders (|z| = {z}). Betroffen sind sechs Modelle.',
 'CHEMBL598527': 'Gehört zum Hauptgerüst der Serie (33 von 39). Die Hexadecylkette ist kein Unterscheidungsmerkmal, 14 Verbindungen der Serie tragen sie. Ein Deskriptor ist bei den übrigen Verbindungen fast konstant, bei dieser Verbindung aber deutlich anders (|z| = {z}). Betroffen ist nur das Modell ohne Auswahl.',
}
KONTROLLE = {'CHEMBL188': 'WIN55212-2, Vergleichssubstanz („higher CB2 selectivity than WIN55212-2"), nicht in der Arbeit hergestellt, Regel Punkt 3'}
SERIEN = [('DS03', 'CHEMBL2016491', 'Garofalo et al. 2012, J Med Chem, (Aryloxy)chinazolin-Harnstoffe, VEGFR-2/EGFR-Assay CHEMBL2020643, IC50, Ersatzregel aktiv ≤ 1 µM'),
          ('DS11', 'CHEMBL1140004', 'Kai et al. 2007, Bioorg Med Chem Lett, CB2-Rezeptor, Assay CHEMBL888169, Ki, Ersatzregel aktiv ≤ 100 nM'),
          ('DS15', 'CHEMBL1159121', 'Lu et al. 2010, Zytotoxizität RH7777, Assay CHEMBL1065834, IC50, Ersatzregel aktiv ≤ 10 µM')]
F = 'Arial'; HF = Font(name=F, bold=True); NF = Font(name=F); TF = Font(name=F, bold=True, size=13)
HEAD = PatternFill('solid', fgColor='D9E1F2'); GRAU = PatternFill('solid', fgColor='D9D9D9'); ORANGE = PatternFill('solid', fgColor='F8CBAD'); GELB = PatternFill('solid', fgColor='FFF2CC')
wb = Workbook(); ws0 = wb.active; ws0.title = 'Legende'
L = [('Numerisch instabile Vorhersagen in DS03, DS11 und DS15 (Hauptanalyse ohne Kontroll- und Referenzsubstanzen)', TF), ('', NF),
     ('Je Serie ein Blatt mit allen Verbindungen, die ChEMBL für den gewählten Assay der Arbeit enthält.', NF),
     ('Grau: Kontroll- oder Referenzsubstanz, nach der Ergänzungsregel ausgeschlossen, nicht in der Rechnung.', NF),
     ('Orange: Ausreißer. Mindestens eine Leave-one-out-Vorhersage eines PLS-Modells über 2 oder unter −1 (Kriterium R1). Begründung in der Spalte „Begründung".', NF),
     ('Gelb: keine instabile Vorhersage, aber labelfreies Kriterium R2 erfüllt (größtes |z| über 100). Diese Verbindungen sind gefährdet, haben aber keine entgleiste Vorhersage erzeugt.', NF), ('', NF),
     ('Vorhersagen: Klasse aktiv = 1, nicht aktiv = 0. Unauffällige Vorhersagen liegen zwischen etwa −0,5 und 1,5.', NF),
     ('Größtes |z|: Abstand der herausgelassenen Verbindung vom Mittelwert der übrigen Verbindungen in deren Standardabweichungen, über alle bereinigten Deskriptoren, ohne Deskriptoren, die bei den übrigen Verbindungen konstant sind (SD < 1e-6). Hängt nicht von den Klassen ab.', NF),
     ('Ursache der Entgleisung: Autoskalierung eines Deskriptors, der bei den Trainingsverbindungen fast konstant ist. Der Abstand der herausgelassenen Verbindung wird dadurch extrem groß und die lineare PLS-Vorhersage schießt über.', NF),
     ('Generisches Gerüst: Bemis-Murcko-Gerüst mit allen Atomen als Kohlenstoff und allen Bindungen als Einfachbindung, wie bei der Auswahl. Spalte „Gerüst in Serie": wie viele Verbindungen der Serie dieses Gerüst teilen.', NF), ('', NF),
     ('Quelle: ChEMBL 37, Zusatzspezifikation vom 24.10.2026 17:20, Ergebnisse der Hauptanalyse (Hauptanalyse_ENDSTAND_mit_Auswertung 25.10.2026 1142).', NF)]
for i, (t, f) in enumerate(L, 1):
    c = ws0.cell(i, 1, t); c.font = f; c.alignment = Alignment(wrap_text=True)
ws0.column_dimensions['A'].width = 150
for dsn, dok, titel in SERIEN:
    info = json.load(open(f'datasets_k/{dsn}_info.json'))['info']
    yk = pd.read_csv(f'datasets_k/{dsn}_y.csv'); Xk = pd.read_csv(f'datasets_k/{dsn}_X.csv')
    idx = {c: i for i, c in enumerate(yk.chembl_id)}
    A = Xk.values.astype(float); cols = list(Xk.columns)
    # label-free max |z| per compound with descriptor
    mz = {}
    for i in range(len(A)):
        tr = np.delete(A, i, 0); mu = tr.mean(0); sd = tr.std(0); ok = sd >= 1e-6
        z = np.abs((A[i, ok] - mu[ok]) / sd[ok]); j = int(np.argmax(z))
        mz[yk.chembl_id[i]] = (float(z[j]), np.array(cols)[ok][j], float(sd[ok][j]))
    # all compounds of the assay from ChEMBL (same query as datasets.py)
    aid = ds.q("select assay_id from assays where chembl_id=?", (info['assay'],))[0][0]
    rows = ds.q("""select molregno, standard_relation, standard_value, standard_units from activities where assay_id=? and src_id='1'
                   and data_validity_comment is null and potential_duplicate='0' and standard_value is not null and molregno is not null""", (str(aid),))
    per = defaultdict(list)
    for m, rel, val, unit in rows: per[ds.parent_of(int(m))].append((rel, float(val), unit))
    keys = {r[0]: (r[1], r[2]) for r in CR.get(dok, [])}
    allc = []
    for m, lst in sorted(per.items()):
        mol, cid = ds.mol_of(m)
        smi = Chem.MolToSmiles(mol) if mol is not None else ''
        werte = ', '.join(f'{r}{v:g}' for r, v, _ in lst)
        allc.append((cid, smi, werte, lst[0][2], mol))
    sc = {}
    for cid, smi, _, _, mol in allc:
        if mol is not None and cid in idx:
            sc[cid] = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(mol)))
    scc = Counter(sc.values())
    ha_med = float(np.median([Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in yk.smiles]))
    ws = wb.create_sheet(dsn)
    ws.cell(1, 1, f'{dsn} – {titel}').font = TF
    ws.cell(2, 1, f"In der Rechnung: {len(yk)} Verbindungen ({int(yk.y.sum())} aktiv, {int((yk.y == 0).sum())} nicht aktiv). Median der Schweratome {ha_med:g}.").font = NF
    head = ['Nr.', 'ChEMBL-ID', 'Bezeichnung im Artikel', 'Name (ChEMBL)', 'SMILES', 'Messwert(e)', 'Einheit', 'Klasse', 'Status', 'Schweratome', 'Gerüst in Serie',
            'Größtes |z| (labelfrei)', 'Deskriptor mit größtem |z|', 'SD dieses Deskriptors bei den übrigen'] + [MOD[m] for m in MOD] + ['Begründung']
    for j, h in enumerate(head, 1):
        c = ws.cell(4, j, h); c.font = HF; c.fill = HEAD; c.alignment = Alignment(wrap_text=True, vertical='top')
    r = 5; nr = 0
    for cid, smi, werte, unit, mol in allc:
        nr += 1; key, name = keys.get(cid, ('', ''))
        drin = cid in idx
        if cid in (info.get('entfernt') or []): status = 'ausgeschlossen: ' + KONTROLLE.get(cid, 'Kontroll- oder Referenzsubstanz')
        elif drin: status = 'in der Rechnung'
        elif mol is None: status = 'nicht verwendet: Struktur nicht verwendbar'
        else: status = 'nicht verwendet: zensierter Wert nicht einstufbar'
        klasse = ('aktiv' if yk.y[idx[cid]] == 1 else 'nicht aktiv') if drin else ''
        preds = [V[dsn]['pred'][m][idx[cid]] if drin else None for m in MOD]
        ausr = drin and any(p > 2 or p < -1 for p in preds)
        z, dname, dsd = mz.get(cid, (None, '', None))
        grund = GRUND.get(cid, '').replace('{z}', f'{z:,.0f}'.replace(',', '.')) if ausr else ''
        if drin and not ausr and z is not None and z > 100:
            grund = f'Kein Ausreißer (alle Vorhersagen unauffällig), aber |z| = {z:,.0f} über 100 (Kriterium R2).'.replace(',', '.')
        vals = [nr, cid, key, name, smi, werte, unit, klasse, status, mol.GetNumHeavyAtoms() if mol is not None else '',
                scc.get(sc.get(cid), '') if drin else '', round(z, 1) if z is not None else '', dname, dsd if dsd is not None else ''] + \
               [round(p, 3) if p is not None else '' for p in preds] + [grund]
        for j, v in enumerate(vals, 1):
            c = ws.cell(r, j, v); c.font = NF; c.alignment = Alignment(wrap_text=(j in (5, 9, len(vals))), vertical='top')
            if not drin and cid in (info.get('entfernt') or []): c.fill = GRAU
            elif ausr: c.fill = ORANGE
            elif drin and z is not None and z > 100: c.fill = GELB
        if ausr:
            for j, p in enumerate(preds, 15):
                if p is not None and (p > 2 or p < -1): ws.cell(r, j).font = Font(name=F, bold=True, color='C00000')
        r += 1
    for j, w in enumerate([5, 15, 16, 22, 45, 16, 9, 10, 34, 9, 9, 11, 22, 12] + [10] * 8 + [70], 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = 'C5'
    print(dsn, 'compounds', nr, 'in computation', len(yk))
wb.save(OUT); print('saved', OUT)
