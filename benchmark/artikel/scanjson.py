import json,re,sys,glob
KEY=re.compile(r'(defined as|cut-?off|threshold|considered (to be )?(active|inactive|potent)|regarded as|classified as|were (deemed|taken) )',re.I)
MICS=re.compile(r'(MIC|IC ?50|IC50|EC50)\s*(values?)?\s*[≤<>=]\s*\d',re.I)
REV=re.compile(r'(this (review|perspective)|in this review|we review|herein we (summar|review)|perspective (provides|presents)|minireview|literature (data|values) (were|was) )',re.I)
SER=re.compile(r'(analog(ue)?s?|derivatives|series of|synthesi[sz]ed)',re.I)
for f in sys.argv[1:]:
    d=json.load(open(f))
    for a in d['articles']:
        t=a.get('full_text') or ''
        print('='*25, a['identifiers']['pmcid'], a['identifiers'].get('doi'))
        print('TITEL:', a['title'][:160])
        print('full-text length:', len(t))
        ab=a.get('abstract') or ''
        print('ABSTRACT:', ab[:450].replace('\n',' '))
        body=re.sub(r'\s+',' ',t)
        sents=re.split(r'(?<=[.]) ',body)
        n=0
        for s in sents:
            if 20<len(s)<350 and (KEY.search(s) or MICS.search(s)):
                print('  *',s[:280]); n+=1
                if n>=12: break
        for s in sents:
            if REV.search(s): print('  TYP>',s[:200]); break
        print('  series indications:', len(SER.findall(body)))
