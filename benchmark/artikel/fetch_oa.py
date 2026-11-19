import urllib.request, gzip, io, re, sys, os
PAIRS={9308496:'PMC9300001_PMC9310000',6675461:'PMC6670001_PMC6679997',6482410:'PMC6480001_PMC6489993',
 4849414:'PMC4840012_PMC4849971',2843756:'PMC2840001_PMC2849939',3298089:'PMC3290003_PMC3300000',
 1855511:'PMC1850758_PMC1859988',3176295:'PMC3170002_PMC3179986'}
for pid,arch in PAIRS.items():
    out=f'PMC{pid}.xml'
    if os.path.exists(out): continue
    url=f'https://ftp.ebi.ac.uk/pub/databases/pmc/oa/{arch}.xml.gz'
    try:
        raw=urllib.request.urlopen(url,timeout=600).read()
    except Exception as e:
        print(pid,'download error',e); continue
    txt=gzip.decompress(raw).decode('utf-8','replace')
    arts=re.split(r'(?=<article[ >])',txt)
    found=[a for a in arts if f'>PMC{pid}<' in a or f'>{pid}<' in a and 'pmcid' in a]
    hit=None
    for a in arts:
        if re.search(rf'<article-id[^>]*pmc[^>]*>\s*(PMC)?{pid}\s*<',a): hit=a;break
    print(pid, len(raw)//10**6,'MB', 'found' if hit else 'NOT found')
    if hit: open(out,'w').write(hit)
