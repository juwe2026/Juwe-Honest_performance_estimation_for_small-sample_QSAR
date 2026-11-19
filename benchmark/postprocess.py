import json, sys, urllib.request, urllib.parse, time
from collections import defaultdict
c = json.load(open('cand_all.json'))
ART = ['MIC', 'Enzym', 'Rezeptor', 'Zytotox']
for x in c:
    x['grund'] = ''
    if not (30 <= x['n_verbindungen'] <= 60):
        x['grund'] = 'Groesse ausserhalb 30-60' if x['n_verbindungen'] < 100 else 'Groesse >= 100 (nur Elternpool Groessenvariation)'
    elif x['diversitaet'] == 'dazwischen':
        x['grund'] = 'Diversitaet weder kongener noch divers'
c.sort(key=lambda x: x['hash'])
# rank within the activity type among the threshold-independent admissible ones
rank = defaultdict(int)
for x in c:
    if not x['grund']:
        rank[x['art']] += 1; x['rang_in_art'] = rank[x['art']]
    else:
        x['rang_in_art'] = None
# parent pool of the size variation: n >= 100, same criteria except size and diversity
ep=defaultdict(int)
for x in c:
    x['elternpool'] = (x['n_verbindungen'] >= 100 and x['art'] in ('MIC', 'Enzym') and x['diversitaet'] != 'dazwischen' and (x['anteil_kleinere_klasse_vorl'] or 0) >= 0.30)
    if x['elternpool']:
        ep[x['art']] += 1; x['rang_elternpool'] = ep[x['art']]
    else: x['rang_elternpool'] = None
# open-access check for the leading candidates (rank <= 40 per type) and parent pool
need = [x for x in c if (x['rang_in_art'] and x['rang_in_art'] <= 40) or (x['rang_elternpool'] and x['rang_elternpool'] <= 10)]
def epmc(pmid, doi):
    qry = f'EXT_ID:{pmid} AND SRC:MED' if pmid else f'DOI:"{doi}"'
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query': qry, 'format': 'json', 'resultType': 'core'})
    for _ in range(3):
        try:
            r = json.load(urllib.request.urlopen(url, timeout=30))
            res = r.get('resultList', {}).get('result', [])
            if not res: return ('unbekannt', None)
            e = res[0]
            return ('ja' if e.get('isOpenAccess') == 'Y' else ('PMC frei' if e.get('inPMC') == 'Y' else 'nein'), e.get('pmcid'))
        except Exception:
            time.sleep(2)
    return ('Fehler', None)
for i, x in enumerate(need):
    if x.get('pmid') or x.get('doi'):
        x['open_access'], x['pmcid'] = epmc(x.get('pmid'), x.get('doi'))
    else:
        x['open_access'], x['pmcid'] = ('keine ID', None)
json.dump(c, open('cand_ranked.json', 'w'), ensure_ascii=False)
from collections import Counter
print(Counter((x['art'], x['diversitaet']) for x in c if not x['grund']))
print('parent pool', Counter(x['art'] for x in c if x['elternpool']))
print('OA', Counter(x.get('open_access') for x in need))
