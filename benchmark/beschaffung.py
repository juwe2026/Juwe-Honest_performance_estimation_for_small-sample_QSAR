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
from collections import defaultdict
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
import select as _s  # noqa
import importlib.util
spec=importlib.util.spec_from_file_location('sel',_ROOT + '/benchmark/select.py'); sel=importlib.util.module_from_spec(spec); spec.loader.exec_module(sel)
cands=json.load(open(_ROOT + '/benchmark/cand_ranked.json')); verdikte=json.load(open(_ROOT + '/benchmark/verdikte.json'))
belegt,prot,bydoc,docs=sel.run(cands,verdikte)
offen=sel.offene_zellen(belegt)
offene_zellen_je_art_div=defaultdict(list)
for a,dv,b in offen: offene_zellen_je_art_div[(a,dv)].append(b)
# for every open cell: next unread publications in hash order,
# whose assays could give a hit in an open band at all
used_tgt=set()
for v in belegt.values():
    used_tgt.add((v['organismus'] or v['ziel']) if v['art']=='MIC' else (v['ziel'] or v['organismus']))
rows=[]
for (a,dv),bands in sorted(offene_zellen_je_art_div.items()):
    n=0
    for d in docs:
        if d in verdikte: continue
        xs=[x for x in bydoc[d] if x['art']==a and x['diversitaet']==dv]
        if not xs: continue
        treffer=[]
        for i,x in enumerate(bydoc[d],1):
            if x['art']!=a or x['diversitaet']!=dv: continue
            bb,frac=sel.qualifies(x)
            tgt=(x['organismus'] or x['ziel']) if a=='MIC' else (x['ziel'] or x['organismus'])
            if bb in bands and tgt not in used_tgt: treffer.append((i,x,bb))
        if not treffer: continue
        n+=1
        i,x,bb=treffer[0]
        rows.append([f'{a} / {dv}', ' + '.join(bands), n, d, x['journal'], x['jahr'], x['titel'], x.get('doi') or '', x.get('pmid'),
                     x.get('open_access') or 'nicht geprueft', len(bydoc[d]), i, x['assay'],
                     (x['organismus'] or x['ziel']) if a=='MIC' else (x['ziel'] or x['organismus']),
                     x['n_verbindungen'], x['aktiv_vorl'], x['nichtaktiv_vorl'], x['anteil_kleinere_klasse_vorl'], bb])
        if n>=8: break
F='Arial'
wb=Workbook(); ws=wb.active; ws.title='Beschaffungsliste'
cols=['Zelle','offene Baender','Prioritaet','Dokument','Journal','Jahr','Titel','DOI','PMID','Open Access',
      'Assays in der Arbeit','erster passender Assay (Nr.)','Assay','Ziel / Organismus','n','vorl. aktiv','vorl. nicht-aktiv','vorl. Anteil','vorl. Band']
ws.append(cols)
H=Font(name=F,bold=True); N=Font(name=F); hf=PatternFill('solid',fgColor='D9E1F2')
for c in ws[1]: c.font=H; c.fill=hf; c.alignment=Alignment(wrap_text=True,vertical='top')
for r in rows: ws.append(r)
for row in ws.iter_rows(min_row=2):
    for c in row: c.font=N
    cc=row[7]
    if cc.value: cc.hyperlink='https://doi.org/'+str(cc.value); cc.font=Font(name=F,color='0563C1',underline='single')
W={'Titel':55,'Zelle':20,'offene Baender':22,'Ziel / Organismus':26,'DOI':26,'Assays in der Arbeit':12,'erster passender Assay (Nr.)':12}
for i,h in enumerate(cols,1): ws.column_dimensions[get_column_letter(i)].width=W.get(h,13)
ws.freeze_panes='A2'; ws.auto_filter.ref=ws.dimensions
# sheet status
w2=wb.create_sheet('Zellbelegung')
w2.append(['Aktivitaetsart','Diversitaet','Band','Status','Assay','Dokument','Ziel / Organismus','n','aktiv','nicht-aktiv','Anteil','Assay-Rang in der Arbeit'])
for c in w2[1]: c.font=H; c.fill=hf; c.alignment=Alignment(wrap_text=True)
grn=PatternFill('solid',fgColor='C6EFCE')
for a in sel.ARTS:
    for dv in sel.DIVS:
        for b in sel.BANDS:
            v=belegt.get((a,dv,b))
            if v:
                w2.append([a,dv,b,'besetzt',v['assay'],v['dokument'],(v['organismus'] or v['ziel']) if a=='MIC' else (v['ziel'] or v['organismus']),v['n_verbindungen'],v['aktiv_vorl'],v['nichtaktiv_vorl'],v['anteil_kleinere_klasse_vorl'],v['_assay_rang']])
                w2.cell(w2.max_row,4).fill=grn
            else: w2.append([a,dv,b,'offen','','','','','','','',''])
for row in w2.iter_rows(min_row=2):
    for c in row: 
        if not c.font.bold: c.font=N
for i,h in enumerate([c.value for c in w2[1]],1): w2.column_dimensions[get_column_letter(i)].width={'Ziel / Organismus':26,'Assay-Rang in der Arbeit':12}.get(h,14)
w2.freeze_panes='A2'
wb.save((_WORK + '/Beschaffungsliste_Artikel 22.10.2026 1320.xlsx'))
print('rows',len(rows),'| per cell', {k:sum(1 for r in rows if r[0]==f'{k[0]} / {k[1]}') for k in offene_zellen_je_art_div})
