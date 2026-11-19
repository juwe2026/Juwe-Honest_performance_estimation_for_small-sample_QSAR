import re,sys,glob,json
KEY=re.compile(r'(defined as|cut-?off|threshold|considered (to be )?(as )?(active|inactive|potent|promising)|regarded as|classified as|deemed|criteri(on|a) for|(were|was) (considered|taken as)|active (compounds?|agents?) (were|was))',re.I)
NUM=re.compile(r'(MIC|MIC\s*50|MIC50|IC\s*50|IC50)\s*(values?\s*)?(of\s*)?[≤<≥>]\s*\d+',re.I)
REV=re.compile(r'(this review|in this review|we review|minireview|this perspective)',re.I)
cand=json.load(open('cand_ranked.json'))
doc2={}
for x in cand: doc2.setdefault(x['dokument'],x)
for f in sorted(glob.glob('txt/*.txt')):
    d=f.split('/')[-1][:-4]; t=open(f,encoding='utf-8',errors='replace').read()
    t=re.sub(r'-\n','',t); body=re.sub(r'\s+',' ',t)
    x=doc2[d]
    print('='*22, d, x.get('rang_in_zelle'), x.get('doi'))
    print('TITELZEILE:', body[:0])
    print(' KOPF:', body[:260])
    hits=[]
    for s in re.split(r'(?<=[.]) ',body):
        if 25<len(s)<330 and (KEY.search(s) or NUM.search(s)): hits.append(s)
    for s in hits[:10]: print('  *',s[:300])
    r=[s for s in re.split(r'(?<=[.]) ',body) if REV.search(s)]
    print('  REVIEW note:', r[0][:150] if r else 'none')
