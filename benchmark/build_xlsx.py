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
import json
from collections import Counter, defaultdict
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
F='Arial'
c=json.load(open(_ROOT + '/benchmark/cand_ranked.json'))
ARTS=['MIC','Enzym','Rezeptor','Zytotox']; DIVS=['kongener','divers']
def band(s):
    if s is None: return 'n/a'
    return '<25 %' if s<0.25 else ('25-35 %' if s<=0.35 else ('35-40 %' if s<0.40 else '40-50 %'))
wb=Workbook()
H=Font(name=F,bold=True); N=Font(name=F); hdrfill=PatternFill('solid',fgColor='D9E1F2'); yel=PatternFill('solid',fgColor='FFFF00')
def sheet(name, header, rows, widths=None, fill_cols=()):
    ws=wb.create_sheet(name)
    ws.append(header)
    for cell in ws[1]: cell.font=H; cell.fill=hdrfill; cell.alignment=Alignment(wrap_text=True,vertical='top')
    for r in rows: ws.append(r)
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.font=N
            if cell.column in fill_cols: cell.fill=yel
    ws.freeze_panes='A2'
    for i,h in enumerate(header,1):
        ws.column_dimensions[get_column_letter(i)].width=(widths or {}).get(h,14)
    ws.auto_filter.ref=ws.dimensions
    return ws
# README
ws=wb.active; ws.title='README'
lines=[('Kandidatenliste ChEMBL-Benchmark, Revision 2',H),
('Stand 21.10.2026, 22:25. Datenquelle ChEMBL 37 (Release 01.05.2026), MySQL-Export, gefiltert eingelesen.',N),
('Grundlage: Vorab-Spezifikation, Entwurf 3 (21.10.2026, 19:30). Es wurden KEINE Modelle, Filter oder Permutationen gerechnet.',N),
('',N),
('Angewandte Kriterien (schwellenunabhaengig):',H),
('- Literaturdaten (Dokumenttyp PUBLICATION, Quelle 1), kein data_validity_comment, potential_duplicate = 0.',N),
('- Aktivitaetsart: MIC gegen Organismus (ohne Parasiten, Viren, Saeuger), Enzym-IC50 (Einzelprotein, confidence >= 8, Proteinklasse Enzym), Rezeptor-Ki (Einzelprotein, confidence >= 8, Membranrezeptor, Ionenkanal oder Kernrezeptor), Zytotoxizitaet IC50/GI50/CC50 gegen Zelllinie.',N),
('- Eine Einheit je Assay. Elternstruktur, einlesbar, keine Gemische, keine Metalle, >= 6 Schweratome. Mehrfachwerte: Median der Logarithmen.',N),
('- Groesse 30 bis 60 Verbindungen (natuerliche Serien). Elternpool Groessenvariation: >= 100 Verbindungen.',N),
('- Diversitaet: kongener = groesstes generisches Murcko-Geruest >= 50 %. Divers = < 25 % und mittlere Tanimoto-Aehnlichkeit (Morgan r=2, 2048 Bit) < 0,30. Dazwischen = ausgeschlossen.',N),
('- Reihenfolge: SHA-256 von "QSAR-R2-" + Dokument-ID + Assay-ID, aufsteigend.',N),
('',N),
('VORLAEUFIG (vorl.): Klassen und Anteile sind mit den Ersatzschwellen gerechnet (MIC <= 8 ug/mL, Enzym <= 1 uM, Rezeptor <= 100 nM, Zytotox <= 10 uM).',H),
('Nach der Spezifikation gilt vorrangig die Schwelle der Autoren. Die endgueltige Klasseneinteilung und das Band stehen erst nach dem Lesen des Artikels fest.',N),
('',N),
('Ablauf:',H),
('1. Artikel der Leseliste in Rangfolge je Zelle lesen und die Schwelle der Autoren eintragen (gelbe Spalten). Freie Artikel liest Claude, die anderen beschaffen die Autoren.',N),
('2. Klassenanteil mit dieser Schwelle berechnen. Erfuellt der Kandidat 25 % und 10 Verbindungen in der kleineren Klasse, erhaelt er seine Zelle (Band nach Anteil), sofern sie noch frei ist.',N),
('3. Weitere Bedingungen: hoechstens ein Datensatz je Dokument und je Ziel/Organismus/Zelllinie.',N),
('4. Sind alle 16 Zellen besetzt, endet die Pruefung. Reicht Stufe 1 nicht, folgen die naechsten Raenge (Blatt Zellen_Top15).',N),
('',N),
('Blaetter: Uebersicht, Leseliste_Stufe1, Zellen_Top15, Elternpool, Alle_Kandidaten.',N)]
for i,(t,f) in enumerate(lines,1):
    ws.cell(i,1,t).font=f
ws.column_dimensions['A'].width=150
# overview
el=[x for x in c if not x['grund']]
rows=[]
for a in ARTS:
    for d in DIVS:
        xs=[x for x in el if x['art']==a and x['diversitaet']==d]
        b=Counter(band(x['anteil_kleinere_klasse_vorl']) for x in xs)
        rows.append([a,d,len(xs),b['<25 %'],b['25-35 %'],b['35-40 %'],b['40-50 %'],b['n/a']])
ws=sheet('Uebersicht',['Aktivitaetsart','Diversitaet','zulaessige Kandidaten','vorl. <25 %','vorl. 25-35 %','vorl. 35-40 %','vorl. 40-50 %','n/a'],rows)
r=ws.max_row+2
ws.cell(r,1,'Ausschlussgruende (schwellenunabhaengig)').font=H
for k,v in Counter(x['grund'] or 'zulaessig' for x in c).most_common():
    r+=1; ws.cell(r,1,k).font=N; ws.cell(r,3,v).font=N
r+=2; ws.cell(r,1,'Kandidaten nach Strukturpruefung gesamt').font=H; ws.cell(r,3,len(c)).font=N
# Leseliste
cols=['Zelle','Rang in Zelle','Assay','Dokument','Journal','Jahr','Titel','DOI','PMID','Open Access','PMCID','Ziel / Organismus / Zelllinie','Typ','n','vorl. aktiv','vorl. nicht-aktiv','vorl. unentschieden','vorl. Anteil kleinere Klasse','vorl. Band','groesstes Geruest','mittl. Tanimoto','Hinweis','Artikel vorhanden (ja/nein)','Schwelle der Autoren (woertlich)','Fundstelle (Seite/Tabelle)','Schwelle angewandt']
W={'Titel':50,'Ziel / Organismus / Zelllinie':28,'Schwelle der Autoren (woertlich)':40,'Hinweis':26,'DOI':28,'Zelle':18}
def row(x,zelle,rang,hint=''):
    doi=x.get('doi') or ''
    return [zelle,rang,x['assay'],x['dokument'],x['journal'],x['jahr'],x['titel'],doi,x.get('pmid'),x.get('open_access'),x.get('pmcid'),
            x['ziel'] if x['art']!='MIC' else (x['organismus'] or x['ziel']),x['standard_type'],x['n_verbindungen'],x['aktiv_vorl'],x['nichtaktiv_vorl'],
            x['n_unentschieden_vorl'],x['anteil_kleinere_klasse_vorl'],band(x['anteil_kleinere_klasse_vorl']),x['groesstes_geruest_anteil'],x['mittl_tanimoto'],hint,'','','','']
les=[]; top=[]
for a in ARTS:
    for d in DIVS:
        xs=[x for x in el if x['art']==a and x['diversitaet']==d][:15]
        seen=set()
        for i,x in enumerate(xs,1):
            hint='gleiches Dokument wie Rang davor' if x['dokument'] in seen else ''
            seen.add(x['dokument'])
            rr=row(x,f'{a} / {d}',i,hint); top.append(rr)
            if i<=6: les.append(rr)
epx=[x for x in c if x.get('rang_elternpool')]
for a in ['MIC','Enzym']:
    for x in [y for y in epx if y['art']==a][:3]:
        les.append(row(x,f'Elternpool {a}',x['rang_elternpool'],'Groessenvariation'))
ncol=len(cols)
ws=sheet('Leseliste_Stufe1',cols,les,W,fill_cols=(ncol-3,ncol-2,ncol-1,ncol))
for rr in ws.iter_rows(min_row=2):
    cell=rr[7]
    if cell.value: cell.hyperlink='https://doi.org/'+cell.value; cell.font=Font(name=F,color='0563C1',underline='single')
sheet('Zellen_Top15',cols[:-4]+['Hinweis'] if False else cols,top,W,fill_cols=(ncol-3,ncol-2,ncol-1,ncol))
sheet('Elternpool',cols,[row(x,f'Elternpool {x["art"]}',x['rang_elternpool'],'') for x in epx],W)
allc=['Hash','Aktivitaetsart','Typ','Einheit','Assay','Dokument','Journal','Jahr','DOI','PMID','Ziel','Organismus','n','vorl. aktiv','vorl. nicht-aktiv','vorl. unentschieden','vorl. Anteil','Geruest-Anzahl','groesstes Geruest','mittl. Tanimoto','Diversitaet','Rang in Art','Ausschlussgrund','Assay-Beschreibung']
sheet('Alle_Kandidaten',allc,[[x['hash'],x['art'],x['standard_type'],x['einheit'],x['assay'],x['dokument'],x['journal'],x['jahr'],x.get('doi'),x.get('pmid'),x['ziel'],x['organismus'],x['n_verbindungen'],x['aktiv_vorl'],x['nichtaktiv_vorl'],x['n_unentschieden_vorl'],x['anteil_kleinere_klasse_vorl'],x['gerueste'],x['groesstes_geruest_anteil'],x['mittl_tanimoto'],x['diversitaet'],x.get('rang_in_art'),x['grund'],x['beschreibung']] for x in c],{'Hash':20,'Assay-Beschreibung':60})
wb.save((_WORK + '/Kandidatenliste_ChEMBL37 21.10.2026 2225.xlsx'))
print(len(les), Counter(r[9] for r in les))
print('to be procured', sum(1 for r in les if r[9] not in ('ja','PMC free')))
