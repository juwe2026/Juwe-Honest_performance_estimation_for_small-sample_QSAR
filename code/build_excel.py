import pandas as pd, numpy as np
from scipy.stats import t as tdist

out=pd.read_pickle('desc_table.pkl')
META=['Monoterpenoid','SMILES']
ACT=[c for c in out.columns if c.startswith(('H. I','P. A','S. A','S. Pneu','S. pyo'))]
DESC=[c for c in out.columns if c not in META+ACT]

REZEPT=[
 ('MolWt','Konstitution','RDKit Descriptors.MolWt','','Mittlere Molmasse (inkl. impliziter H)'),
 ('C_Count','Konstitution','Atomzählung Z=6','','Anzahl Kohlenstoffatome'),
 ('O_Count','Konstitution','Atomzählung Z=8','','Anzahl Sauerstoffatome'),
 ('HeavyAtomCount','Konstitution','mol.GetNumHeavyAtoms()','','Nicht-H-Atome'),
 ('RingCount','Konstitution','CalcNumRings','','SSSR-Ringzahl'),
 ('NumAromaticRings','Struktur','CalcNumAromaticRings','','RDKit-Default-Aromatizität'),
 ('NumRotatableBonds','Struktur','CalcNumRotatableBonds (Strict)','','Strikte Definition'),
 ('FractionCSP3','Struktur','CalcFractionCSP3','','Anteil sp3-C'),
 ('CeqC_DoubleBonds','Struktur','SMARTS-Zählung','[CX3]=[CX3]','Aliphatische C=C'),
 ('NumStereoCenters','Struktur','FindMolChiralCenters(includeUnassigned=True)','','Zugewiesene + mögliche Stereozentren'),
 ('NumAliphaticRings','Struktur','CalcNumAliphaticRings','','Nicht-aromatische Ringe'),
 ('DegreeUnsaturation','Struktur','(2C+2+N-H-X)/2','','Aus Summenformel (H nach AddHs)'),
 ('NumHDonors','Funktion','CalcNumHBD','','RDKit-verfeinerte HBD'),
 ('NumHAcceptors','Funktion','CalcNumHBA','','RDKit-verfeinerte HBA'),
 ('Aldehyde_flag','Funktion','SMARTS (0/1)','[CX3H1](=O)[#6]','Aldehyd'),
 ('Ketone_flag','Funktion','SMARTS (0/1)','[#6][CX3](=[OX1])[#6]','Keton (inkl. Chinon)'),
 ('CarboxAcid_flag','Funktion','SMARTS (0/1)','[CX3](=O)[OX2H1]','Carbonsäure'),
 ('Phenol_flag','Funktion','SMARTS (0/1)','[c][OX2H1]','Aromatisches OH (inkl. Tropolon)'),
 ('Ether_flag','Funktion','SMARTS (0/1)','[OD2]([#6;!$([CX3]=[OX1])])[#6;!$([CX3]=[OX1])]','Ether; Ester ausgeschlossen'),
 ('EsterLacton_flag','Funktion','SMARTS (0/1)','[CX3](=[OX1])[OX2][#6]','Ester + Lacton'),
 ('Alcohol_flag','Funktion','SMARTS (0/1)','[CX4][OX2H]','Aliphatisches OH'),
 ('MolLogP','Physikochemie','Crippen.MolLogP','','Wildman-Crippen, neutrale Form'),
 ('TPSA','Physikochemie','Descriptors.TPSA','','Ertl; nur N,O'),
 ('MolMR','Physikochemie','Crippen.MolMR','','Molare Refraktivität'),
 ('LabuteASA','Physikochemie','Descriptors.LabuteASA','','2D-Oberflächennäherung'),
]
REZ=pd.DataFrame(REZEPT,columns=['Deskriptor','Gruppe','Methode/RDKit','SMARTS','Konvention/Notiz'])

X=out[DESC].astype(float); n=len(X)
P=X.corr('pearson'); S=X.corr('spearman')
def pmat(C):
    C=C.values.copy(); np.fill_diagonal(C,0)
    with np.errstate(divide='ignore',invalid='ignore'): tt=C*np.sqrt((n-2)/(1-C**2))
    pv=2*tdist.sf(np.abs(tt),n-2); pv[np.abs(C)>=1]=0
    return pv
maxabs=np.maximum(P.abs().values,S.abs().values)
adj=((P.abs().values>0.7)&(pmat(P)<0.05))|((S.abs().values>0.7)&(pmat(S)<0.05))
np.fill_diagonal(adj,False)

# Seltene Binär-Flags (<4 Treffer) -> feste Singletons (instabile Korrelationen)
is01=[set(np.unique(X[d]))<= {0,1} for d in DESC]
rare=[i for i,d in enumerate(DESC) if is01[i] and X[d].sum()<4]
for i in rare: adj[i,:]=False; adj[:,i]=False
print('Feste Singletons (seltene Flags <4):',[DESC[i] for i in rare])

# Greedy Cluster; Repräsentant = höchste Zentralität (mittlere |r| im Cluster), unsupervised
remaining=set(range(len(DESC))); groups=[]
partners={i:set(np.where(adj[i])[0]) for i in range(len(DESC))}
while remaining:
    cnt={i:len(partners[i]&remaining) for i in remaining}; mx=max(cnt.values())
    if mx==0: break
    hub=sorted([i for i in remaining if cnt[i]==mx],key=lambda i:(-cnt[i],i))[0]
    mem=[hub]+sorted(partners[hub]&remaining)
    cent={i:np.mean([maxabs[i,j] for j in mem if j!=i]) for i in mem}
    rep=sorted(mem,key=lambda i:(-cent[i], is01[i], i))[0]  # höchste Zentralität, Nicht-Flag bevorzugt
    groups.append((rep,mem,cent)); remaining-=set(mem)
singletons=sorted(remaining)
reduced=[DESC[rep] for rep,_,_ in groups]+[DESC[i] for i in singletons]
print('Cluster:',len(groups),'| Singletons:',len(singletons),'| reduziert:',len(reduced))
for rep,mem,cent in groups:
    print(f'  [{DESC[rep]}] <- {[DESC[i] for i in mem if i!=rep]}')
print('  Singletons:',[DESC[i] for i in singletons])

red_rows=[]
for rep,mem,cent in groups:
    others=[DESC[i] for i in sorted([m for m in mem if m!=rep],key=lambda i:-cent[i])]
    red_rows.append({'Repräsentant (behalten)':DESC[rep],'Typ':'Cluster',
                     'redundante Partner (|r|>0.7)':', '.join(others),'Clustergröße':len(mem)})
for i in singletons:
    red_rows.append({'Repräsentant (behalten)':DESC[i],
                     'Typ':'Singleton'+(' (seltenes Flag)' if i in rare else ''),
                     'redundante Partner (|r|>0.7)':'','Clustergröße':1})
RED=pd.DataFrame(red_rows)

fn='/mnt/user-data/outputs/20260706_31MT_LitDeskriptoren_Analyse.xlsx'
with pd.ExcelWriter(fn,engine='openpyxl') as w:
    out[META+ACT+DESC].to_excel(w,sheet_name='Deskriptoren',index=False)
    REZ.to_excel(w,sheet_name='Rechenvorschriften',index=False)
    P.round(4).to_excel(w,sheet_name='Pearson_Corr')
    S.round(4).to_excel(w,sheet_name='Spearman_Corr')
    RED.to_excel(w,sheet_name='Reduktion',index=False)
np.save('/home/claude/reduced_desc.npy',np.array(reduced,dtype=object))
out[META+ACT+DESC].to_pickle('/home/claude/desc_full.pkl')
print('\nExcel-Basis geschrieben. Reduzierter Satz ('+str(len(reduced))+'):',reduced)
