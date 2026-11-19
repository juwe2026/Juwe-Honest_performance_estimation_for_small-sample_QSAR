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
import json, importlib.util, re
from collections import Counter
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
spec=importlib.util.spec_from_file_location('sel',_ROOT + '/benchmark/select.py'); sel=importlib.util.module_from_spec(spec); spec.loader.exec_module(sel)
c=json.load(open(_ROOT + '/benchmark/cand_ranked.json')); v=json.load(open(_ROOT + '/benchmark/verdikte.json'))
belegt,prot,bydoc,docs=sel.run(c,v)
F='Arial'; H=Font(name=F,bold=True); N=Font(name=F)
hf=PatternFill('solid',fgColor='D9E1F2'); grn=PatternFill('solid',fgColor='C6EFCE'); red=PatternFill('solid',fgColor='FFC7CE')
wb=Workbook()
def sh(name,cols,rows,W=None):
    ws=wb.create_sheet(name); ws.append(cols)
    for cc in ws[1]: cc.font=H; cc.fill=hf; cc.alignment=Alignment(wrap_text=True,vertical='top')
    for r in rows: ws.append(r)
    for row in ws.iter_rows(min_row=2):
        for cc in row: cc.font=N
    for i,h in enumerate(cols,1): ws.column_dimensions[get_column_letter(i)].width=(W or {}).get(h,14)
    ws.freeze_panes='A2'; ws.auto_filter.ref=ws.dimensions; return ws
# README
ws=wb.active; ws.title='README'
L=[('Versuchsplan ChEMBL-Benchmark – endgueltige Belegung',H),
('Stand 22.10.2026, 17:05. Datenquelle ChEMBL 37 (Release 01.05.2026). Grundlage: Vorab-Spezifikation, Entwurf 9.',N),
('Es wurden KEINE Modelle, Filter oder Permutationen gerechnet. Alle Zahlen stammen aus der Auswahl.',N),('',N),
('Alle 16 Zellen des Versuchsplans sind besetzt: 4 Aktivitaetsarten x 2 Diversitaetsstufen x 2 Baender.',H),
('Je Publikation genau ein Datensatz. 15 Datensaetze wurden im ersten Durchgang vergeben, 1 im zweiten Durchgang mit Zielwiederholung (Mycobacterium tuberculosis).',N),('',N),
('Klassendefinition: KEINE der 48 gelesenen Arbeiten nennt eine Schwelle fuer aktiv gegenueber nicht-aktiv.',H),
('Es gilt daher durchgehend die Ersatzregel: MIC <= 8 ug/mL, Enzym-IC50 <= 1 uM, Rezeptor-Ki <= 100 nM, Zytotoxizitaet <= 10 uM.',N),('',N),
('Ausgeschlossene Publikationen (Blatt Ausschluesse): Uebersichtsarbeiten, Arbeiten mit zusammengetragenen Werten, Arbeiten mit nach Aktivitaet vorselektierten Verbindungen.',N),('',N),
('Blaetter: Zellbelegung, Datensaetze, Ausschluesse, Auswahlprotokoll, Alle_Kandidaten.',N)]
for i,(t,f) in enumerate(L,1): ws.cell(i,1,t).font=f
ws.column_dimensions['A'].width=150
def tgt(x): return (x['organismus'] or x['ziel']) if x['art']=='MIC' else (x['ziel'] or x['organismus'])
# cell assignment
rows=[]
for a in sel.ARTS:
    for dv in sel.DIVS:
        for b in sel.BANDS:
            x=belegt.get((a,dv,b))
            rows.append([a,dv,b,x['dokument'],x.get('doi'),x['assay'],tgt(x),x['n_verbindungen'],x['aktiv_vorl'],x['nichtaktiv_vorl'],
                         x['anteil_kleinere_klasse_vorl'],x['groesstes_geruest_anteil'],x['mittl_tanimoto'],x['_assay_rang'],x.get('_durchgang',1),x['journal'],x['jahr']])
ws=sh('Zellbelegung',['Aktivitaetsart','Diversitaet','Band','Dokument','DOI','Assay','Ziel / Organismus / Zelllinie','n','aktiv','nicht-aktiv','Anteil kleinere Klasse','groesstes Geruest','mittl. Tanimoto','Assay-Nr. in der Arbeit','Durchgang','Journal','Jahr'],rows,
   {'Ziel / Organismus / Zelllinie':28,'DOI':30,'Anteil kleinere Klasse':12,'Assay-Nr. in der Arbeit':12,'Journal':20})
for row in ws.iter_rows(min_row=2): row[3].fill=grn
# data sets DS01..DS16
rows=[]
for i,((a,dv,b),x) in enumerate(sorted(belegt.items()),1):
    rows.append([f'DS{i:02d}',a,dv,b,x['dokument'],x.get('doi'),x['assay'],tgt(x),x['standard_type'],x['einheit'],x['n_verbindungen'],
                 x['aktiv_vorl'],x['nichtaktiv_vorl'],x['anteil_kleinere_klasse_vorl'],x['gerueste'],x['groesstes_geruest_anteil'],x['mittl_tanimoto'],
                 'Ersatzregel (Autoren nennen keine Schwelle)','Originalarbeit mit eigener Serie'])
sh('Datensaetze',['ID','Aktivitaetsart','Diversitaet','Band','Dokument','DOI','Assay','Ziel','Typ','Einheit','n','aktiv','nicht-aktiv','Anteil','Gerueste','groesstes Geruest','mittl. Tanimoto','Schwelle','Art der Arbeit'],rows,
   {'DOI':30,'Ziel':26,'Art der Arbeit':28,'Schwelle':34})
# exclusions
TXT={'review':'Uebersichtsarbeit','keine_serie':'Werte aus anderen Publikationen zusammengetragen','vorselektiert':'Verbindungen nach Aktivitaet vorselektiert (Screening-Treffer)','methoden':'Methoden-/Assay-Arbeit ohne eigene Serie','unlesbar':'Volltext nicht verfuegbar'}
rows=[]
for d,info in sorted(v.items()):
    if info['art']=='original': continue
    xs=bydoc.get(d) or []
    j=xs[0] if xs else None
    rows.append([d,info.get('doi'),TXT.get(info['art'],info['art']),(j['journal'] if j else ''),(j['jahr'] if j else ''),(j['titel'][:120] if j else '')])
ws=sh('Ausschluesse',['Dokument','DOI','Grund','Journal','Jahr','Titel'],rows,{'DOI':30,'Grund':46,'Titel':60})
for row in ws.iter_rows(min_row=2): row[2].fill=red
# selection log
rows=[[d,a or '',g] for d,a,g in prot]
sh('Auswahlprotokoll',['Dokument','Assay','Entscheidung'],rows,{'Entscheidung':60})
# all candidates (compact)
sh('Alle_Kandidaten',['Hash','Aktivitaetsart','Assay','Dokument','DOI','Ziel','Organismus','n','aktiv','nicht-aktiv','Anteil','Diversitaet','Ausschlussgrund'],
   [[x['hash'][:16],x['art'],x['assay'],x['dokument'],x.get('doi'),x['ziel'],x['organismus'],x['n_verbindungen'],x['aktiv_vorl'],x['nichtaktiv_vorl'],x['anteil_kleinere_klasse_vorl'],x['diversitaet'],x['grund']] for x in c],
   {'DOI':28,'Ziel':24,'Ausschlussgrund':30})
wb.save((_WORK + '/Versuchsplan_ChEMBL_Zellbelegung 22.10.2026 1705.xlsx'))
print('cells',len(belegt),'| exclusions',len(rows) if False else sum(1 for i in v.values() if i['art']!='original'),'| log rows',len(prot))
print(Counter(x['diversitaet'] for x in belegt.values()), Counter(k[2] for k in belegt))
