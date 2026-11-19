"""Raw MIC replicates of the 31 monoterpenoids from the laboratory workbook raw.xlsx (QSAR_RAW_MIC_XLSX, not part of the
archive), written to raw_mic_31.json in the output folder. The replicate values are given in Additional file 2 of the
article, sheet MIC_replicates_31."""
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
import openpyxl, re, numpy as np, json, warnings, sys
sys.path.insert(0,(_ROOT + '/code')); import core
warnings.filterwarnings('ignore')
ws=openpyxl.load_workbook(_osp.environ.get('QSAR_RAW_MIC_XLSX', _WORK + '/raw.xlsx'),data_only=True)['Original']
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
json.dump(OUT,open((_OUT + '/raw_mic_31.json'),'w'),ensure_ascii=False,indent=0)
