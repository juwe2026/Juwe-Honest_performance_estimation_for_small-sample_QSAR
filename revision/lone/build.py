import csv, re, json
from rdkit import Chem
from rdkit.Chem import Descriptors
# D-amino acid backbone: N[C@H](R)C(=O)...  (L would be [C@@H] in this writing order)
SIDE={'k':'CCCCN','a':'C','f':'Cc1ccccc1','s':'CO','t':'[C@@H](C)O','y':'Cc1ccc(O)cc1','v':'C(C)C','r':'CCCNC(=N)N',
      '1-nal':'Cc1cccc2ccccc12','nle':'CCCC','dap':'CN'}
def parse(seq):
    s=seq[:-4] if seq.endswith('-NH2') else seq
    res=[];i=0
    while i<len(s):
        if s[i]=='(':
            j=s.index(')',i); res.append(s[i+1:j]); i=j+1
        else:
            res.append(s[i]); i+=1
    return res
def smiles(seq):
    res=parse(seq); out=[]
    for k,r in enumerate(res):
        sd=SIDE[r]
        if r=='t': out.append('N[C@@H]([C@@H](C)O)C(=O)') if False else out.append('N[C@H]([C@H](C)O)C(=O)')
        else: out.append('N[C@H](%s)C(=O)'%sd)
    return ''.join(out)+'N'   # C-terminal amide
rows=[]
with open('t1.csv',encoding='utf-8-sig') as f:
    rd=list(csv.reader(f,delimiter=';'))
hdr=[i for i,r in enumerate(rd) if r and r[0]=='No'][0]
cols=rd[hdr]
for r in rd[hdr+1:]:
    if not r or not r[0] or r[0]=='Strain key': break
    rows.append(dict(zip(cols,r)))
print(len(rows))
for r in rows:
    sm=smiles(r['Sequence']); m=Chem.MolFromSmiles(sm)
    assert m is not None, r
    r['SMILES']=Chem.MolToSmiles(m); r['MW']=round(Descriptors.MolWt(m),2)
    r['nres']=len(parse(r['Sequence']))
json.dump(rows,open('rows.json','w'),indent=1)
for r in rows[:3]+rows[-2:]: print(r['No'],r['Sequence'],r['nres'],r['MW'])
