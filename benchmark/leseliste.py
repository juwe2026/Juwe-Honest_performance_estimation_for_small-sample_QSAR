import json, urllib.request, urllib.parse, time
from collections import Counter
c=json.load(open('cand_ranked.json'))
def epmc(pmid, doi):
    qry = f'EXT_ID:{pmid} AND SRC:MED' if pmid else f'DOI:"{doi}"'
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query': qry, 'format': 'json', 'resultType': 'core'})
    for _ in range(3):
        try:
            r = json.load(urllib.request.urlopen(url, timeout=30)); res = r.get('resultList', {}).get('result', [])
            if not res: return ('unbekannt', None)
            e = res[0]
            return ('ja' if e.get('isOpenAccess') == 'Y' else ('PMC frei' if e.get('inPMC') == 'Y' else 'nein'), e.get('pmcid'))
        except Exception: time.sleep(2)
    return ('Fehler', None)
el=[x for x in c if not x['grund']]
for art in ['MIC','Enzym','Rezeptor','Zytotox']:
    for div in ['kongener','divers']:
        xs=[x for x in el if x['art']==art and x['diversitaet']==div]
        for i,x in enumerate(xs[:15],1):
            x['rang_in_zelle']=i
            if x.get('open_access') in (None,'Fehler'):
                x['open_access'],x['pmcid']=epmc(x.get('pmid'),x.get('doi'))
ep=[x for x in c if x.get('rang_elternpool')]
for x in ep:
    if x.get('open_access') in (None,'Fehler'): x['open_access'],x['pmcid']=epmc(x.get('pmid'),x.get('doi'))
json.dump(c,open('cand_ranked.json','w'),ensure_ascii=False)
print(Counter(x.get('open_access') for x in c if x.get('rang_in_zelle')))
