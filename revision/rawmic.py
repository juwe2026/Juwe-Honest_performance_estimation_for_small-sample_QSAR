import openpyxl, re, numpy as np, json, warnings, sys
sys.path.insert(0,'QSAR_selection_leakage/code'); import core
warnings.filterwarnings('ignore')
ws=openpyxl.load_workbook('raw.xlsx',data_only=True)['Original']
rows=list(ws.iter_rows(values_only=True)); hdr=rows[0]
PAT={'HI':r'H\.I\.','PA':r'P\.A\.','SA':r'S\. aureus','SPneu':r'S\.? ?Pneu','SPyo':r'S\.? ?Pyo'}
def cols(k,phase): return [j for j,h in enumerate(hdr) if h and phase in str(h) and re.search(PAT[k],str(h))]
def norm(s): return re.sub(r'[^a-z0-9]','',s.lower().replace('α','a').replace('β','b'))
rawidx={norm(r[1]):i for i,r in enumerate(rows) if i>0 and r[1]}
def val(x):
    if x is None or (isinstance(x,str) and not x.strip()): return None
    if isinstance(x,str):
        s=x.strip().replace('˃','>')
        if s.startswith('>'): return ('cens',float(s[1:]))
        return ('val',float(s))
    return ('val',float(x))
OUT={}
for c in core.COMPOUNDS:
    n=norm(c); key=[k for k in rawidx if k==n or k.startswith(n[:8])]; assert len(key)==1,(c,key)
    r=rows[rawidx[key[0]]]; OUT[c]={'raw_name':r[1]}
    for k in PAT:
        for ph in ['broth','agar']:
            vs=[val(r[j]) for j in cols(k,ph)]; vs=[v for v in vs if v]
            # coded: censored >1024 -> 2048 ; numeric 2048 also means >1024 (S. aureus sheet)
            num=np.array([2048.0 if (t=='cens' or x>=2048) else x for t,x in vs])
            OUT[c][k+'_'+ph]=dict(n=len(num),n_inhib=int((num<=1024).sum()),mean=float(num.mean()) if len(num) else None,
                                 median=float(np.median(num)) if len(num) else None,values=num.tolist())
json.dump(OUT,open('out/raw_mic_31.json','w'),ensure_ascii=False,indent=0)
