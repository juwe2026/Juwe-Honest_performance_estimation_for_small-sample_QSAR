"""Excel working workbook: per data set all analysed compounds with the control and reference substances marked.
Basis: supplementary specification DRAFT 2 of 24.10.2026, 06:00. Literature from the ChEMBL API (document), keys from compound_record.
The cell texts of this working workbook are in German (internal working document). The decisions are reported in English
in Additional file 12 (sheet Excluded_compounds) and Additional file 8.
"""
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
import json, sqlite3, subprocess, re, html
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
def gscaf(m):
    try: return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))
    except Exception: return ''
PRUEF = {  # structurally related: check by EW (24.10.2026, 08:45)
 'DS09': {'CHEMBL2443262': 'Oliceridin (3): Serie aus TRV130 abgeleitet, aber kein gemeinsames Grundgerüst',
          'CHEMBL5415557': 'TRV734 (5): Analogon von Oliceridin, kein gemeinsames Grundgerüst'},
 'DS10': {'CHEMBL841': 'Loperamid (3): Vorlage der MOR-Hälfte der Hybride, gemeinsames Grundgerüst mit Serienmitgliedern',
          'CHEMBL1201294': 'Diphenoxylat (4): Vorlage der MOR-Hälfte der Hybride, gemeinsames Grundgerüst mit Serienmitgliedern',
          'CHEMBL651': 'Methadon (6): Vorlage der MOR-Hälfte der Hybride, gemeinsames Grundgerüst mit Serienmitgliedern'},
 'DS12': {'CHEMBL56': '8-OH-DPAT: gleiche pharmakologische Klasse (5-HT1A-Agonist), anderes Grundgerüst als die Serie'},
}
ORANGE = PatternFill('solid', fgColor='FCE4D6')

BASIS = (_ROOT + '/benchmark')
CR = json.load(open(_ROOT + '/benchmark/kontrollen/compound_records.json'))
META = json.load(open(_ROOT + '/benchmark/kontrollen/docs_meta.json'))
META['CHEMBL4190169'].update(authors='Ragab FA, Heiba HI, El-Gazzar MG, Abou-Seri SM, El-Sabbagh WA, El-Hazek RM.',
                              title='Synthesis of novel thiadiazole derivatives as selective COX-2 inhibitors',
                              doi=META['CHEMBL4190169'].get('doi') or '10.1039/c6md00367b')
db = sqlite3.connect(f'file:{BASIS}/chembl37_subset.db?mode=ro', uri=True)   # read-only: a missing database is an error, not a new empty file

AUS = {  # chembl_id: (point of the rule, evidence in the article)
 'DS02': {'CHEMBL118': ('3', '„with the reference drugs celecoxib and diclofenac"'),
          'CHEMBL139': ('3', '„with the reference drugs celecoxib and diclofenac"')},
 'DS05': {'CHEMBL64': ('3', '„anti-tuberculosis drug isoniazid (MIC 0.01 μg/mL), which was used as a control"'),
          'CHEMBL44884': ('3', '„comparable to the drug ethambutol (0.8 μg/mL)", nicht in der Arbeit hergestellt')},
 'DS06': {c: ('3', 'Antibiotikum zur Charakterisierung der Resistenzmuster der Stämme (Table 2, Table S2), gehört nicht zur Serie der PNA-Peptid-Konjugate')
          for c in ['CHEMBL8', 'CHEMBL127', 'CHEMBL702', 'CHEMBL148', 'CHEMBL158', 'CHEMBL44354', 'CHEMBL174', 'CHEMBL177',
                    'CHEMBL186', 'CHEMBL376140', 'CHEMBL1440', 'CHEMBL1730']},
 'DS07': {'CHEMBL64': ('3', '„Isoniazid (INH), which is a first line antituberculosis drug, was used as a standard", INH in der MIC-Tabelle')},
 'DS08': {'CHEMBL106': ('3', '„fluconazole as reference drug"')},
 'DS09': {'CHEMBL2443262': ('3', '„Reference compounds 3, 5, and 6"'), 'CHEMBL5415557': ('3', '„Reference compounds 3, 5, and 6"')},
 'DS10': {c: ('3', '„Reference compounds, including fentanyl, 3, 4, 5, and 6, … are reported in the table for useful comparisons"')
          for c in ['CHEMBL54', 'CHEMBL596', 'CHEMBL651', 'CHEMBL841', 'CHEMBL1201294']},
 'DS11': {'CHEMBL188': ('3', 'Vergleichssubstanz („higher CB2 selectivity than WIN55212-2"), nicht in der Arbeit hergestellt')},
 'DS12': {'CHEMBL56': ('3', 'Vergleichssubstanz („5 times more potent than 8-OH-DPAT"), nicht in der Arbeit hergestellt')},
 'DS13': {'CHEMBL53463': ('3', '„for the positive control doxorubicin and the precursors lapachone"')},
 'DS14': {'CHEMBL185': ('3 und 4', 'Vergleichssubstanz, Wert mit hochgestelltem Literaturverweis (8)')},
 'DS16': {'CHEMBL53463': ('3', '„compared with the positive control doxorubicin"')},
 'DS21': {'CHEMBL91': ('3', '„systemic reference drugs such as … miconazole …"'), 'CHEMBL157101': ('3', '„systemic reference drugs such as … ketoconazole …"')},
 'DS22': {},
 'DS23': {'CHEMBL5179655': ('4', '„Data for Leucettines L90−L92 are from ref 65"')},
 'DS24': {},
}
GRENZ = {
 'DS01': {'CHEMBL4297507': 'Delgocitinib (ent-60): in dieser Arbeit entdeckt und hergestellt, Endprodukt der Serie (Punkt 1)'},
 'DS05': {'CHEMBL214920': 'Isoxyl (43): in der Arbeit synthetisiert („42 and 43 (isoxyl)"), in den SAR-Tabellen der Serie, bekannter Arzneistoff (Punkt 1, Grenzfall)'},
 'DS10': {'CHEMBL346331': 'Normethadon (19): „we prepared the desmethyl-methadone 19" (Punkt 1)'},
 'DS12': {'CHEMBL349003': 'Leitverbindung 1 (U93385), Titel: „Analogs of … (1, U93385)" (Punkt 2)'},
 'DS13': {c: 'Lapachon: „the precursors lapachone", Ausgangsstoff der Derivate (Punkt 2)' for c in ['CHEMBL15192', 'CHEMBL441441', 'CHEMBL394882']},
 'DS15': {},
}


def schluessel(doc, cid):
    ks = [k for m, k, n in CR[doc] if m == cid]
    if ks:
        return ' / '.join(ks)
    m = db.execute('select molregno from molecule_dictionary where chembl_id=?', (cid,)).fetchone()
    if not m:
        return ''
    out = []
    for (k,) in db.execute('select molregno from molecule_hierarchy where parent_molregno=?', (m[0],)).fetchall():
        kid = db.execute('select chembl_id from molecule_dictionary where molregno=?', (k,)).fetchone()[0]
        ks = [kk for mm, kk, nn in CR[doc] if mm == kid]
        if not ks:
            r = json.loads(subprocess.run(['curl', '-s', '--max-time', '30',
                                           f'https://www.ebi.ac.uk/chembl/api/data/compound_record.json?molecule_chembl_id={kid}'],
                                          capture_output=True, text=True).stdout)
            ks = [x['compound_key'] for x in r['compound_records'] if x['document_chembl_id'] == doc]
        out += [f'{x} (Salzform {kid})' for x in ks]
    return ' / '.join(out)


def zitat(doc):
    m = META[doc]
    t = re.sub(r'<[^>]+>', '', html.unescape(m['title'] or '')).rstrip('.')
    seiten = f"{m['first_page']}" + (f"–{m['last_page']}" if m.get('last_page') and m['last_page'] != m['first_page'] else '')
    return f"{(m['authors'] or '').rstrip('.')}. {t}. {m['journal']} {m['year']};{m['volume']}:{seiten}."


ARIAL = 'Arial'
F = lambda **k: Font(name=ARIAL, **k)
ROT = PatternFill('solid', fgColor='F4CCCC')
GELB = PatternFill('solid', fgColor='FFF2CC')
GRAU = PatternFill('solid', fgColor='D9D9D9')
GRUEN = PatternFill('solid', fgColor='D9EAD3')
dünn = Side(style='thin', color='BFBFBF')
RAND = Border(left=dünn, right=dünn, top=dünn, bottom=dünn)

wb = Workbook()
ue = wb.active; ue.title = 'Übersicht'
dss = [f'DS{i:02d}' for i in range(1, 17)] + ['DS21', 'DS22', 'DS23', 'DS24']
zeilen_ue = []
for ds in dss:
    info = json.load(open(f'{BASIS}/datasets/{ds}_info.json'))['info']
    doc = info['dokument']
    t = pd.read_csv(f'{BASIS}/datasets/{ds}_y.csv')
    ws = wb.create_sheet(ds)
    aus = AUS.get(ds, {}); gr = GRENZ.get(ds, {})
    if ds == 'DS15':
        for j, c in enumerate(t.chembl_id):
            k = schluessel(doc, c)
            if k in ('1a', '1b'):
                gr[c] = f'{k}: in der Arbeit hergestellt, zugleich „Compounds 1a and 1b were chosen as control compounds" (Punkt 1, Grenzfall)'
    if ds in ('DS22', 'DS24'):
        for j, c in enumerate(t.chembl_id):
            k = schluessel(doc, c)
            if k == 'Amphotericin B':
                aus[c] = ('3', '„systemic reference drugs such as … amphotericin B"')
            if k == 'L90':
                aus[c] = ('4', '„Data for Leucettines L90−L92 are from ref 65"')
    teil = ds in ('DS21', 'DS22', 'DS23', 'DS24')
    ws['A1'] = f'{ds} – {" / ".join(info["zelle"])}'; ws['A1'].font = F(bold=True, size=13)
    kopf = [('Literatur', zitat(doc)), ('DOI', f"https://doi.org/{META[doc]['doi']}"), ('ChEMBL-Dokument', doc),
            ('Assay / Ziel', f"{info['assay']} / {info['ziel']}"), ('Einheit der Messwerte', info['einheit'])]
    if teil:
        kopf.append(('Hinweis', f"Teilmenge (n = {info['n']}) aus dem Elterndatensatz {info['assay']}. Geht nicht in H1 bis H4 ein. Ob Teilmengen bereinigt werden, entscheidet EW."))
    ms = [Chem.MolFromSmiles(x) for x in t.smiles]; fps = [GEN.GetFingerprint(m) for m in ms]; scs = [gscaf(m) for m in ms]
    kern = [j for j, c in enumerate(t.chembl_id) if c not in aus]
    simmax, gleich = [], []
    for j in range(len(t)):
        andere = [k for k in kern if k != j]
        simmax.append(max(DataStructs.BulkTanimotoSimilarity(fps[j], [fps[k] for k in andere])))
        gleich.append(sum(scs[k] == scs[j] for k in andere))
    nn = [simmax[j] for j in kern]
    npr = len([c for c in aus if c in PRUEF.get(ds, {})])
    status = ('Keine Änderung notwendig: keine Kontroll- oder Referenzsubstanz in diesem Datensatz.' if not aus
              else f'{len(aus)} Verbindung(en) auszuschließen (rot markiert)' + (f', davon {npr} zur Prüfung (orange, strukturell verwandt).' if npr else '.'))
    kopf.append(('Ähnlichkeit innerhalb der Serie', f'nächster Nachbar (Morgan, Radius 2, Tanimoto) im Median {np.median(nn):.2f}, Minimum {min(nn):.2f}. Vergleichswert für Spalte K.'))
    for r, (a, b) in enumerate(kopf, start=2):
        ws.cell(r, 1, a).font = F(bold=True); ws.cell(r, 2, b).font = F()
        ws.cell(r, 2).alignment = Alignment(wrap_text=True, vertical='top')
        ws.merge_cells(start_row=r, start_column=2, end_row=r, end_column=8)
    ws.row_dimensions[2].height = 45
    ws.cell(len(kopf) + 1, 2).alignment = Alignment(wrap_text=True, vertical='top')
    r0 = len(kopf) + 2
    ws.cell(r0, 1, 'Status').font = F(bold=True)
    ws.cell(r0, 2, status).font = F(bold=True, color='9C0006' if aus else '006100')
    ws.cell(r0, 2).fill = ROT if aus else GRUEN
    ws.merge_cells(start_row=r0, start_column=2, end_row=r0, end_column=8)
    hz = r0 + 2
    spalten = ['Nr.', 'ChEMBL-ID', 'Schlüssel im Artikel', 'ChEMBL-Name', 'Relation', 'Messwert', 'Klasse',
               'Entscheidung', 'Punkt der Regel', 'Beleg / Begründung', 'max. Tanimoto zur Serie',
               'gleiches Grundgerüst wie … Serienmitglieder']
    for c, h in enumerate(spalten, start=1):
        z = ws.cell(hz, c, h); z.font = F(bold=True); z.fill = GRAU; z.border = RAND
        z.alignment = Alignment(wrap_text=True, vertical='center')
    for j, row in t.iterrows():
        c = row['chembl_id']
        pn = db.execute('select pref_name from molecule_dictionary where chembl_id=?', (c,)).fetchone()
        pn = pn[0] if pn and pn[0] else ''
        if c in aus and c in PRUEF.get(ds, {}):
            ent, pkt, beleg, fill = 'ausschließen (zur Prüfung)', aus[c][0], aus[c][1] + ' | Prüfhinweis: ' + PRUEF[ds][c], ORANGE
        elif c in aus:
            ent, pkt, beleg, fill = 'ausschließen', aus[c][0], aus[c][1], ROT
        elif c in gr:
            ent, pkt, beleg, fill = 'bleibt (Grenzfall)', '', gr[c], GELB
        else:
            ent, pkt, beleg, fill = 'bleibt', '', '', None
        werte = [j + 1, c, schluessel(doc, c), pn, row['relation'], float(row['wert']),
                 'aktiv' if row['y'] == 1 else 'nicht aktiv', ent, pkt, beleg, round(simmax[j], 2), int(gleich[j])]
        for k, v in enumerate(werte, start=1):
            z = ws.cell(hz + 1 + j, k, v); z.font = F(); z.border = RAND
            z.alignment = Alignment(wrap_text=(k == 10), vertical='top')
            if k in (11, 12): z.alignment = Alignment(horizontal='center', vertical='top')
            if fill: z.fill = fill
            if k == 5:
                z.value = str(v); z.data_type = 's'
        ws.cell(hz + 1 + j, 6).number_format = '0.###'
    n = len(t); a, b = hz + 1, hz + n
    s = b + 2
    rechn = [
        ('Verbindungen bisher', f'=COUNTA(B{a}:B{b})'),
        ('davon auszuschließen', f'=COUNTIF(H{a}:H{b},"ausschließen*")'),
        ('Verbindungen nachher', f'=B{s}-B{s+1}'),
        ('aktiv nachher', f'=COUNTIF(G{a}:G{b},"aktiv")-COUNTIFS(G{a}:G{b},"aktiv",H{a}:H{b},"ausschließen*")'),
        ('nicht aktiv nachher', f'=COUNTIF(G{a}:G{b},"nicht aktiv")-COUNTIFS(G{a}:G{b},"nicht aktiv",H{a}:H{b},"ausschließen*")'),
        ('kleinere Klasse nachher', f'=MIN(B{s+3},B{s+4})'),
        ('Anteil kleinere Klasse', f'=IF(B{s+2}>0,B{s+5}/B{s+2},0)'),
        ('zulässig (n ≥ 30, kleinere Klasse ≥ 10 und ≥ 25 %)' if not teil else 'zulässig (für Teilmengen nicht geprüft)',
         f'=IF(AND(B{s+2}>=30,B{s+5}>=10,B{s+6}>=0.25),"ja","nein")' if not teil else 'entfällt'),
    ]
    for i, (lab, f) in enumerate(rechn):
        ws.cell(s + i, 1, lab).font = F(bold=True); ws.cell(s + i, 2, f).font = F()
    ws.cell(s + 6, 2).number_format = '0.0%'
    for col, w in zip('ABCDEFGHIJKL', [24, 16, 22, 22, 9, 11, 11, 22, 9, 70, 12, 16]):
        ws.column_dimensions[col].width = w
    ws.cell(s + 9, 1, 'Spalte H darf geändert werden: „ausschließen" oder „bleibt". Die Zählung oben und die Übersicht rechnen sich dann neu.').font = F(italic=True, size=9)
    ws.freeze_panes = ws.cell(hz + 1, 1)
    zeilen_ue.append((ds, info, doc, s, teil, len(aus)))

# overview
ue['A1'] = 'Kontroll- und Referenzsubstanzen im ChEMBL-Benchmark – Übersicht'; ue['A1'].font = F(bold=True, size=13)
ue['A2'] = ('Grundlage: Zusatzspezifikation ENTWURF 2 vom 24.10.2026, 06:00 (noch nicht festgeschrieben). '
            'Die primäre Analyse bleibt unverändert. Rot = auszuschließen, strukturfremd. Orange = auszuschließen nach Entwurf 2, aber strukturell verwandt, zur Prüfung durch EW. '
            'Gelb = Grenzfall, bleibt in der Serie. Spalte K/L je Datensatz: Ähnlichkeit zur übrigen Serie.')
ue['A2'].font = F(italic=True); ue.merge_cells('A2:J2'); ue['A2'].alignment = Alignment(wrap_text=True); ue.row_dimensions[2].height = 30
kopf = ['Datensatz', 'Zelle', 'Literatur', 'bisher n', 'auszuschließen', 'nachher n', 'kleinere Klasse', 'Anteil', 'zulässig', 'Status', 'davon zur Prüfung']
for c, h in enumerate(kopf, start=1):
    z = ue.cell(4, c, h); z.font = F(bold=True); z.fill = GRAU; z.border = RAND
for i, (ds, info, doc, s, teil, na) in enumerate(zeilen_ue):
    r = 5 + i
    m = META[doc]
    werte = [ds, ' / '.join(info['zelle']), f"{(m['authors'] or '').split(',')[0]} et al., {m['journal']} {m['year']}",
             f"='{ds}'!B{s}", f"='{ds}'!B{s+1}", f"='{ds}'!B{s+2}", f"='{ds}'!B{s+5}", f"='{ds}'!B{s+6}", f"='{ds}'!B{s+7}",
             'keine Änderung' if na == 0 else 'Änderung']
    for c, v in enumerate(werte, start=1):
        z = ue.cell(r, c, v); z.font = F(color='008000') if isinstance(v, str) and v.startswith('=') else F(); z.border = RAND
        if c == 10: z.fill = GRUEN if na == 0 else ROT
    ue.cell(r, 8).number_format = '0.0%'
    z = ue.cell(r, 11, f'=COUNTIF(\'{ds}\'!H:H,"ausschließen (zur Prüfung)")'); z.font = F(color='008000'); z.border = RAND
r = 5 + len(zeilen_ue) + 1
ue.cell(r, 1, 'Summe auszuschließen, natürliche Serien').font = F(bold=True)
ue.cell(r, 5, '=SUM(E5:E20)').font = F(bold=True)
ue.cell(r + 1, 1, 'Serien mit Änderung').font = F(bold=True)
ue.cell(r + 1, 5, '=COUNTIF(J5:J20,"Änderung")').font = F(bold=True)
ue.cell(r + 3, 1, 'Regel (Entwurf 2): Zur Serie gehören die neu vorgestellten Verbindungen und ihre Ausgangs-, Leit- oder Vorläuferstrukturen (Punkte 1 und 2). '
        'Nicht zur Serie gehören fremde Arzneistoffe, klinische Kandidaten oder Werkzeugsubstanzen, die zum Vergleich, als Referenz, Standard oder Positivkontrolle '
        'mitgemessen wurden (Punkt 3), und Werte, die der Artikel aus einer anderen Publikation übernimmt (Punkt 4). Entscheidend ist die Herkunft, nicht das Wort.').font = F()
ue.merge_cells(start_row=r + 3, start_column=1, end_row=r + 3, end_column=10)
ue.cell(r + 3, 1).alignment = Alignment(wrap_text=True, vertical='top'); ue.row_dimensions[r + 3].height = 60
ue.cell(r + 4, 1, 'Quellen: Literaturangaben aus der ChEMBL-API (document), DS02 ergänzt aus der ersten PDF-Seite. Schlüssel aus ChEMBL compound_record. Belege aus den Verlags-PDFs im Ordner „Artikel".').font = F(italic=True, size=9)
ue.merge_cells(start_row=r + 4, start_column=1, end_row=r + 4, end_column=10)
for col, w in zip('ABCDEFGHIJK', [12, 30, 34, 9, 14, 10, 14, 9, 9, 16, 14]):
    ue.column_dimensions[col].width = w
out = (_WORK + '/Kontrollsubstanzen_je_Datensatz 24.10.2026 0900.xlsx')
wb.save(out); print(out)
