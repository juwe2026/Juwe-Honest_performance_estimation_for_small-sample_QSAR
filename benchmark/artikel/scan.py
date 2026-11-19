import sys, re, urllib.request, time, json, html
KEY = re.compile(r'(active|inactive|activ\w+ (was|were) defined|cut-?off|threshold|considered\s+\w*\s*active|MIC\s*[<≤>≥]|IC50\s*[<≤>≥]|regarded as|classified as|potent\w*\s+(was|were)?\s*defined|criteri\w+)', re.I)
TYPE = re.compile(r'(this (review|perspective)|in this review|we review|herein we summar|perspective (provides|presents)|minireview)', re.I)
def get(p):
    for _ in range(4):
        try:
            x = urllib.request.urlopen(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{p}/fullTextXML', timeout=60).read().decode('utf-8','replace')
            if len(x) > 2000: return x
        except Exception as e:
            pass
        time.sleep(3)
    return None
for p in sys.argv[1:]:
    x = get(p)
    print('='*20, p, 'XML' if x else 'NO FULL TEXT')
    if not x: continue
    art = re.search(r'article-type="([^"]+)"', x)
    print('article-type:', art.group(1) if art else '?')
    body = re.sub(r'<[^>]+>', ' ', x)
    body = html.unescape(re.sub(r'\s+', ' ', body))
    sents = re.split(r'(?<=[.;]) ', body)
    seen = set()
    for s in sents:
        if len(s) > 400 or len(s) < 25: continue
        if KEY.search(s) and s not in seen:
            seen.add(s); print(' *', s[:300])
        if len(seen) > 25: break
    for s in sents:
        if TYPE.search(s): print(' TYP>', s[:250]); break
