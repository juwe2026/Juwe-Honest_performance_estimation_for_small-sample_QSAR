import numpy as np, pandas as pd, openpyxl
from openpyxl.styles import Font, PatternFill
nlv1=np.load('lit_nlv1.npy',allow_pickle=True).item()
cv=pd.read_pickle('lit_results.pkl'); two=pd.read_pickle('lit_twostage.pkl')
fn='/mnt/user-data/outputs/20260706_31MT_LitDeskriptoren_Analyse.xlsx'

# 2208er-Ergebnisse (aus früheren Sitzungen)
biased={'H.I.':0.43,'S.A.':0.50,'S.Pneu':0.56,'S.Pyo':0.52,'P.A.':0.88}
nested2208={'H.I.':(-0.14,0.15),'S.A.':(0.16,0.013),'S.Pneu':(0.19,0.010),'S.Pyo':(-0.11,0.13),'P.A.':(0.40,0.0015)}
order=['H.I.','S.A.','S.Pneu','S.Pyo','P.A.']; full={'H.I.':'H. influenzae','S.A.':'Staph. aureus','S.Pneu':'S. pneumoniae','S.Pyo':'S. pyogenes','P.A.':'P. aeruginosa'}
cvm={r['Erreger']:r for _,r in cv.iterrows()}

# PLS-DA_reduziert neu: nLV=1 (ehrlich) + nLV-CV (Obergrenze)
rows=[]
for s in order:
    Q1,A1,pQ1,pA1=nlv1[s]; c=cvm[s]
    rows.append({'Erreger':full[s],'aktiv':c['aktiv'],
                 'Q2 (nLV=1, EHRLICH)':round(Q1,3),'AUC (nLV=1)':round(A1,3),'p(Q2) nLV=1':round(pQ1,4),'p(AUC) nLV=1':round(pA1,4),
                 'nLV (CV)':c['nLV'],'Q2 (nLV-CV, Obergrenze)':c['Q2'],'AUC (nLV-CV)':c['AUC'],
                 'Konfusion nLV-CV (TP/FN/FP/TN)':c['Konfusion']})
RES=pd.DataFrame(rows)

# Drei-Wege-Vergleich
cmp=[]
for s in order:
    Q1,A1,pQ1,pA1=nlv1[s]; nq,npv=nested2208[s]
    note=''
    if s=='P.A.': note='Artefakt (3 Aktive) - ausschliessen'
    elif s=='S.Pyo': note='schwach/nicht signifikant in allen Ansaetzen'
    else: note='Literatur-Set holt ehrliches, signifikantes Signal zurueck'
    cmp.append({'Erreger':full[s],
                '(1) 2208 verzerrt Q2':biased[s],
                '(2) 2208 genestet Q2':nq,'(2) p(Q2)':npv,
                '(3) Literatur-13 Q2 (nLV=1, ehrlich)':round(Q1,3),'(3) p(Q2)':round(pQ1,4),
                'Interpretation':note})
CMP=pd.DataFrame(cmp)

with pd.ExcelWriter(fn,engine='openpyxl',mode='a',if_sheet_exists='replace') as w:
    RES.to_excel(w,sheet_name='PLS-DA_reduziert',index=False)
    CMP.to_excel(w,sheet_name='Vergleich_3Ansaetze',index=False)

# Fazit-Blatt mit Kommentaren
wb=openpyxl.load_workbook(fn)
if 'Fazit' in wb.sheetnames: del wb['Fazit']
ws=wb.create_sheet('Fazit'); H=Font(bold=True,size=12); fill=PatternFill('solid',fgColor='DDEEDD'); r=1
def put(t,sec=False):
    global r; c=ws.cell(row=r,column=1,value=t)
    if sec: c.font=H;c.fill=fill
    r+=1
put('FAZIT – Literaturbasierte Deskriptorauswahl (leakage-frei)',sec=True)
for t in [
 '1) Kernbefund: Mit 13 unsupervised-reduzierten, literaturbasierten Deskriptoren zeigen H. influenzae,',
 '   S. aureus und S. pneumoniae ein EHRLICH signifikantes, moderates Signal (Q2=0.29-0.39, p<0.01)',
 '   - selbst bei nLV=1. Genau diese Endpunkte hatte der 2208-Deskriptor-Ansatz genestet fast (H.I.,',
 '   S.Pyo) bzw. teils (S.A., S.Pneu) verloren.',
 '2) Diagnostischer Wert (wie vorab erwartet): Der Ansatz TRENNT die Ursachen. Das Versagen der',
 '   2208-Pipeline war nicht (nur) fehlendes Signal, sondern Overfitting/Leakage der supervidierten',
 '   Selektion aus zu vielen Kandidaten bei n=31. Ein kleines, theoriebasiertes, leakage-freies Set',
 '   macht das reale (moderate) Struktur-Wirkungs-Signal sichtbar.',
 '3) Warum hier eine EINFACHE Permutation ehrlich ist: Die Reduktion nutzt nur Deskriptor-Deskriptor-',
 '   Korrelation (kein Y) -> leakage-frei. Es gibt keinen supervidierten Selektionsschritt zu nesten.',
 '   apparent ~ kreuzvalidiert (kein grosser Einbruch) - im Gegensatz zum 2208-Ansatz.',
 '4) S. pyogenes bleibt schwach (Q2=0.06, n.s.) trotz der meisten Aktiven (11) -> endpunkt-spezifisch,',
 '   nicht von der Aktivenzahl getrieben. Biologisch interpretierbar.',
 '5) P. aeruginosa (3 Aktive) bleibt ein Overfit-Artefakt (AUC=1.0) und wird nur deskriptiv gefuehrt.',
 '6) Multiplizitaet: Unter Bonferroni (alpha=0.01) bestehen H.I. (0.001), S.Pneu (0.0025), S.A. (p(Q2)',
 '   =0.003; p(AUC)=0.018 grenzwertig). In der Publikation angeben.',
 '',
 'ZUR VIP-FRAGE (PLS-DA mit vielen Deskriptoren, dann nur hohe VIP behalten, neu rechnen):',
 '- Als Selektionsmethode fuer die GUETESCHAETZUNG nur legitim, wenn die VIP-Selektion IN jeder',
 '  CV-Fold wiederholt wird (genestet). Sonst ist es dieselbe Leakage-Falle wie beim Effektgroessen-',
 '  Filter - nur mit VIP als Kriterium.',
 '- Beleg aus Blatt VIP_zweistufig: naive VIP-Selektion (auf allen Daten) > genestete VIP-Selektion',
 '  bei Q2 durchgaengig (z.B. H.I. 0.29 -> 0.16; S.Pneu 0.53 -> 0.31). Der Unterschied IST der',
 '  Leakage-Betrag.',
 '- Empfehlung: VIP zur INTERPRETATION nutzen (welche Deskriptoren tragen), nicht als zweite Selektion',
 '  zur Aufhuebschung der berichteten Guete. Mit einem guten a-priori-Set ist die zweite Stufe ohnehin',
 '  meist unnoetig - der reduzierte Satz modelliert bereits ehrlich.',
 '',
 'HINWEIS nLV: nLV=1 ist die konservative, ehrliche Basis. Die nLV-CV-Werte (nLV 3-4) sind hoeher,',
 'enthalten aber einen Optimismus durch die Komponentenwahl per Q2; als Obergrenze zu lesen.',
]:
    put(t)
ws.column_dimensions['A'].width=100
wb.save(fn)
print('Final geschrieben. Blätter:', openpyxl.load_workbook(fn).sheetnames)
