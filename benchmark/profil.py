import re,sys,glob,os
KEY=re.compile(r'(defined as (active|inactive|a hit)|cut-?off (value|of|for)|considered (to be )?(as )?(active|inactive|potent|a hit)|regarded as (active|inactive)|classified as (active|inactive)|criteri(on|a) for (activity|a hit)|(compounds?|analogues?|derivatives?) (with|showing) (MIC|IC50|Ki|EC50)[^.]{0,40}(were|was) (considered|regarded|classified|defined))',re.I)
SER=re.compile(r'(a series of|we (designed|synthesi[sz]ed|report)|were (designed and )?synthesi[sz]ed|novel .{0,30}(derivatives|analogues|analogs))',re.I)
def find(doi):
    pat=doi.replace('/','_').replace('+','+')
    for f in glob.glob('txt/*.txt'):
        b=os.path.basename(f)[:-4]
        n=b.replace(':','_')
        if n.lower()==pat.lower(): return f
        if doi.startswith('10.1021/acs.jmedchem.'):
            tail=doi.split('acs.jmedchem.')[1]
            if b.lower() in ('jm'+tail.lower(),): return f
    return None
for doi in sys.argv[1:]:
    f=find(doi)
    print('='*20, doi, '->', os.path.basename(f) if f else 'FILE NOT FOUND')
    if not f: continue
    t=re.sub(r'-\n','',open(f,encoding='utf-8',errors='replace').read())
    b=re.sub(r'\s+',' ',t)
    print(' KOPF:', b[:330])
    hits=[s for s in re.split(r'(?<=[.]) ',b) if 25<len(s)<320 and KEY.search(s)]
    print(' THRESHOLD candidates:', len(hits))
    for s in hits[:5]: print('   *', s[:300])
    ser=[s for s in re.split(r'(?<=[.]) ',b) if 25<len(s)<260 and SER.search(s)]
    print(' SERIE:', ser[0][:230] if ser else 'no indication')
