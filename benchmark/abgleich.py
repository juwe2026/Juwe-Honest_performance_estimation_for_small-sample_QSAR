import json, re, importlib.util
spec=importlib.util.spec_from_file_location('sel','select.py'); sel=importlib.util.module_from_spec(spec); spec.loader.exec_module(sel)
cands=json.load(open('cand_ranked.json')); verdikte=json.load(open('verdikte.json'))
belegt,prot,bydoc,docs=sel.run(cands,verdikte)
offen=sel.offene_zellen(belegt)
from collections import defaultdict
offzell=defaultdict(list)
for a,dv,b in offen: offzell[(a,dv)].append(b)
used_tgt=set()
for v in belegt.values():
    used_tgt.add((v['organismus'] or v['ziel']) if v['art']=='MIC' else (v['ziel'] or v['organismus']))
# regenerate the procurement list (identical logic as beschaffung.py)
liste=defaultdict(list)   # (art,div) -> [(prio, doc, doi, assay, ziel, band)]
for (a,dv),bands in sorted(offzell.items()):
    n=0
    for d in docs:
        if d in verdikte: continue
        treffer=[]
        for i,x in enumerate(bydoc[d],1):
            if x['art']!=a or x['diversitaet']!=dv: continue
            bb,frac=sel.qualifies(x)
            tgt=(x['organismus'] or x['ziel']) if a=='MIC' else (x['ziel'] or x['organismus'])
            if bb in bands and tgt not in used_tgt: treffer.append((i,x,bb))
        if not treffer: continue
        n+=1; i,x,bb=treffer[0]
        liste[(a,dv)].append((n,d,(x.get('doi') or '').lower(),tgt if False else ((x['organismus'] or x['ziel']) if a=='MIC' else (x['ziel'] or x['organismus'])),bb,x['anteil_kleinere_klasse_vorl'],x['journal'],x['jahr']))
        if n>=8: break
# existing files
namen=open('dateinamen.txt').read().split()
def norm(f):
    s=f[:-4] if f.endswith('.pdf') else f
    s=s.replace(':','/')
    if s.startswith('10.'):
        s=re.sub(r'^(10\.\d{4,5})_', r'\1/', s)
    else:
        m=re.fullmatch(r'jm(\d)c(\d{5})', s)
        s='10.1021/acs.jmedchem.'+m.group(1)+'c'+m.group(2) if m else '10.1021/'+s
    return s.lower()
vorhanden={norm(f):f for f in namen}
print('files:',len(vorhanden))
alle_doi={}
for k,v in liste.items():
    for t in v: alle_doi[t[2]]=(k,t)
tref=set(vorhanden)&set(alle_doi)
print('hits in the procurement list:',len(tref),'von',len(alle_doi))
print('\nfiles without assignment to the procurement list:')
for d in sorted(set(vorhanden)-tref): print('  ',d, vorhanden[d])
print()
for k,v in sorted(liste.items()):
    da=[t for t in v if t[2] in vorhanden]; fehlt=[t for t in v if t[2] not in vorhanden]
    print(f"== {k[0]} / {k[1]}  open: {' + '.join(offzell[k])}  available {len(da)} of {len(v)}")
    for t in v:
        mark='OK ' if t[2] in vorhanden else 'FEHLT'
        print(f"   {mark} Prio {t[0]}  {t[2]:38s} {t[4]:11s} {str(t[5]):6s} {t[3][:26]:28s} {t[6]} {t[7]}")
