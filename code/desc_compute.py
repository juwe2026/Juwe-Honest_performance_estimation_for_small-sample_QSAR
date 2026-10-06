import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 6 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import pandas as pd, numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

# SMARTS-Definitionen (getrennt Aldehyd/Keton; Ether ohne Ester)
SM = {
 'Aldehyde_flag'  : '[CX3H1](=O)[#6]',
 'Ketone_flag'    : '[#6][CX3](=[OX1])[#6]',
 'CarboxAcid_flag': '[CX3](=O)[OX2H1]',
 'Phenol_flag'    : '[c][OX2H1]',
 'Ether_flag'     : '[OD2]([#6;!$([CX3]=[OX1])])[#6;!$([CX3]=[OX1])]',  # schließt Ester-O aus
 'EsterLacton_flag': '[CX3](=[OX1])[OX2][#6]',   # Ester + Lacton (O an C gebunden, nicht H)
 'Alcohol_flag'   : '[CX4][OX2H]',               # aliphatisches -OH
 'CeqC_count'     : '[CX3]=[CX3]',   # aliphatische C=C
}
SMOBJ={k:Chem.MolFromSmarts(v) for k,v in SM.items()}

def dou(mol):
    # Degree of Unsaturation = (2C+2+N-H-X)/2  (X=Halogene)
    cnt={}
    mh=Chem.AddHs(mol)
    for a in mh.GetAtoms(): cnt[a.GetSymbol()]=cnt.get(a.GetSymbol(),0)+1
    C=cnt.get('C',0); H=cnt.get('H',0); N=cnt.get('N',0)
    X=sum(cnt.get(x,0) for x in ['F','Cl','Br','I'])
    return (2*C+2+N-H-X)/2

def compute(smiles):
    mol=Chem.MolFromSmiles(smiles)           # Sanitize + RDKit-Default-Aromatizität
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    logp,mr=Crippen.MolLogP(mol),Crippen.MolMR(mol)
    n_atoms=lambda z:sum(1 for a in mol.GetAtoms() if a.GetAtomicNum()==z)
    def nmatch(key): return len(mol.GetSubstructMatches(SMOBJ[key],uniquify=True))
    # Revision (October 2026): CIP-based count (useLegacyImplementation=True, the RDKit default), as in the
    # revised Additional file 3. The submitted file used useLegacyImplementation=False, which also counts the two
    # ring-bridgehead atoms of 1,8-cineole (2 instead of 0); every other value of the 25 descriptors is the same.
    stereo=len(Chem.FindMolChiralCenters(mol,includeUnassigned=True,useLegacyImplementation=True))
    d={
     # Grundlegende Konstitution
     'MolWt'            : Descriptors.MolWt(mol),
     'C_Count'          : n_atoms(6),
     'O_Count'          : n_atoms(8),
     'HeavyAtomCount'   : mol.GetNumHeavyAtoms(),
     'RingCount'        : rdMolDescriptors.CalcNumRings(mol),
     'NumAromaticRings' : rdMolDescriptors.CalcNumAromaticRings(mol),
     'NumRotatableBonds': rdMolDescriptors.CalcNumRotatableBonds(mol,rdMolDescriptors.NumRotatableBondsOptions.Strict),
     # Sättigung / Struktur
     'FractionCSP3'     : rdMolDescriptors.CalcFractionCSP3(mol),
     'CeqC_DoubleBonds' : nmatch('CeqC_count'),
     'NumStereoCenters' : stereo,
     'NumAliphaticRings': rdMolDescriptors.CalcNumAliphaticRings(mol),
     'DegreeUnsaturation': dou(mol),
     # Funktionelle Gruppen
     'NumHDonors'       : rdMolDescriptors.CalcNumHBD(mol),
     'NumHAcceptors'    : rdMolDescriptors.CalcNumHBA(mol),
     'Aldehyde_flag'    : int(nmatch('Aldehyde_flag')>0),
     'Ketone_flag'      : int(nmatch('Ketone_flag')>0),
     'CarboxAcid_flag'  : int(nmatch('CarboxAcid_flag')>0),
     'Phenol_flag'      : int(nmatch('Phenol_flag')>0),
     'Ether_flag'       : int(nmatch('Ether_flag')>0),
     'EsterLacton_flag' : int(nmatch('EsterLacton_flag')>0),
     'Alcohol_flag'     : int(nmatch('Alcohol_flag')>0),
     # Physikochemie
     'MolLogP'          : logp,
     'TPSA'             : Descriptors.TPSA(mol),
     'MolMR'            : mr,
     'LabuteASA'        : Descriptors.LabuteASA(mol),
    }
    return d

if __name__=='__main__':
    xl=pd.ExcelFile(_os.environ.get("QSAR_ROUTEA_XLSX", _os.path.join(_DATA, XLSX)))
    df=pd.read_excel(xl,'31MT_5Bac_AllDes')
    meta=['Monoterpenoid','SMILES']+[c for c in df.columns if c.startswith(('H. I','P. A','S. A','S. Pneu','S. pyo'))]
    rows=[compute(s) for s in df['SMILES']]
    D=pd.DataFrame(rows); order=list(rows[0].keys())
    out=pd.concat([df[meta].reset_index(drop=True),D],axis=1)
    out.to_pickle(_os.path.join(_RES, 'desc_table.pkl'))
    print('Descriptors:',len(order))
    # Heikle Fälle prüfen
    print('\n== Functional-group flags (which molecules) ==')
    for flag in ['Aldehyde_flag','Ketone_flag','CarboxAcid_flag','Phenol_flag','Ether_flag']:
        names=out.loc[out[flag]==1,'Monoterpenoid'].tolist()
        print(f'{flag:16}: {len(names)} -> {names}')
    print('\n== Aromatic rings (quinone check) ==')
    print(out.loc[out['NumAromaticRings']>0,['Monoterpenoid','NumAromaticRings','Phenol_flag']].to_string(index=False))
    print('\n== Suspicious: molecules with O but NO functional-group flag ==')
    noflag=out[(out['O_Count']>0)&(out[['Aldehyde_flag','Ketone_flag','CarboxAcid_flag','Phenol_flag','Ether_flag']].sum(axis=1)==0)]
    print(noflag[['Monoterpenoid','SMILES','O_Count']].to_string(index=False) if len(noflag) else '  none')
