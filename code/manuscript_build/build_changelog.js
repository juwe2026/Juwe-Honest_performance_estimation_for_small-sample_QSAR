const d = require("docx");
const { Document, Packer, Paragraph, TextRun, AlignmentType, HeadingLevel } = d;
const fs = require("fs");
const FONT = "Calibri";
const ch = [];
const push = x => ch.push(x);
function H1(t){ return new Paragraph({ spacing:{before:280,after:140}, children:[new TextRun({text:t,bold:true,size:30,font:FONT,color:"1F4E79"})]}); }
function H2(t){ return new Paragraph({ spacing:{before:220,after:100}, children:[new TextRun({text:t,bold:true,size:24,font:FONT,color:"2E74B5"})]}); }
function P(t){ return new Paragraph({ spacing:{after:120,line:264}, alignment:AlignmentType.JUSTIFIED, children:[new TextRun({text:t,size:21,font:FONT})]}); }
function It(nr,title,text){
  push(new Paragraph({ spacing:{before:120,after:40}, children:[new TextRun({text:nr+"  "+title,bold:true,size:21,font:FONT})]}));
  push(new Paragraph({ spacing:{after:120,line:264}, alignment:AlignmentType.JUSTIFIED, children:[new TextRun({text:text,size:21,font:FONT})]}));
}

push(new Paragraph({ spacing:{after:60}, children:[new TextRun({text:"Manuskript-Aenderungsdatei",bold:true,size:34,font:FONT,color:"1F4E79"})]}));
push(new Paragraph({ spacing:{after:60}, children:[new TextRun({text:"Werle, Rondevaldova und Werle - Journal of Cheminformatics",bold:true,size:24,font:FONT})]}));
push(new Paragraph({ spacing:{after:220}, children:[new TextRun({text:"Fassung: Werle_JCheminform_revised_Claude_2026-08-10.docx (28 Seiten) mit Supplement (22 Seiten). Stand dieser Datei: 11.08.2026, Runde zum genetischen Algorithmus.",italics:true,size:21,font:FONT})]}));
push(P("Diese Runde ergaenzt einen Vergleich mit einer Wrapper-Methode und behebt zwei formale Defekte in der Literaturliste, die dabei aufgefallen sind."));

push(H1("1. Neue Rechnung: genetischer Algorithmus"));
It("1.1","Was gerechnet wurde.","Ein GA nach dem Bauplan der OPERA-Arbeit: binaeres Chromosom je Deskriptor, Population 30, Crossover 0.5, Mutation 0.01, 50 Generationen, Fitness gleich kreuzvalidiertes Q2 mit Strafterm fuer die Modellgroesse. Gesucht wurde auf den 182 label-freien Repraesentanten, nicht im vollen Pool, der rund 10 hoch 665 Teilmengen umfasst. Zwei Staemme: H. influenzae und, auf euren Wunsch, S. aureus als der Endpunkt mit dem besten Q2 unserer Filtermethode.");
It("1.2","Ergebnis.","H. influenzae: Der GA berichtet selbst Q2 = 0.64 bei 26 Deskriptoren; dieselbe Auswahl per LOO bewertet ergibt 0.45; wird die Suche in jeder Fold wiederholt, bleiben -0.59. S. aureus: selbst berichtet 0.72 bei 24 Deskriptoren, LOO 0.70, genestet -0.08. Der Absturz betraegt damit 1.23 und 0.80 gegenueber 0.42 bei unserem Effektgroessenfilter.");
It("1.3","Warum das der Kernaussage nutzt.","Der Wrapper leckt zwei- bis dreimal staerker als der Filter, und am staerksten ausgerechnet bei S. aureus, dem einzigen Endpunkt, fuer den unsere Filtermethode ein reproduzierbares Signal findet. Die Suchfreiheit und nicht die Wahl des Algorithmus bestimmt also, wie weit ein Verfahren die eigene Schaetzung aufblaeht. Der GA waehlte zudem in jeder Fold einen anderen Satz, im Median 23 Deskriptoren bei einer Spanne von 13 bis 33.");
It("1.4","Faire Einordnung der zitierten Arbeit.","Absatz und Legende halten ausdruecklich fest, dass dies keine Kritik an der OPERA-Arbeit ist. Dort wurde vor dem GA ein externer Testsatz von 25 Prozent abgetrennt und die berichtete Guete stammt von diesem; die Leckagefrage stellt sich dort nicht. Der Befund illustriert nur, was passiert, wenn dieselbe Suche ohne externen Testsatz durchgefuehrt wird.");
It("1.5","Umfang der Rechnung und ihre Grenzen.","Ein GA-Lauf dauert 9.2 Sekunden. Das vollstaendige Protokoll der zitierten Arbeit mit 100 Laeufen in zwei Runden waere fuer einen Stamm 31 Minuten, genestet 79 Stunden, mit 2,000 Permutationen rund 18 Jahre. Unsere Zahlen stammen daher aus einem Lauf je Fold und sind indikativ. Das steht so in der Legende von Table S12.");

push(H1("2. Umsetzung im Dokument"));
It("2.1","Neuer Absatz im Diskussionsteil.","Er steht unmittelbar vor dem Absatz zum Verhaeltnis von R2 und Q2 und knuepft an den dort begonnenen Groessenvergleich der Leckagequellen an, den er nach oben verlaengert.");
It("2.2","Neue Table S12 im Supplement.","Sie stellt fuer beide Staemme gegenueber: die vom Verfahren selbst berichtete Guete, dieselbe Auswahl per LOO, die genestete Guete, den Absturz, die Modellgroesse je Fold und zum Vergleich die Werte unserer Route A1.");
It("2.3","Neue Referenz.","Mansouri und Kollegen 2018, OPERA-Modelle, J Cheminform 10:10.");

push(H1("3. Zwei formale Defekte behoben"));
It("3.1","Referenz zum OECD-Leitfaden war unzitiert.","Bei der Verlagerung zweier Abschnitte ins Supplement war die einzige Zitierstelle des OECD-Validierungsleitfadens verloren gegangen. Ein Eintrag im Literaturverzeichnis, der nirgends zitiert wird, ist ein formaler Mangel. Die Conclusions enthalten jetzt einen Satz, der die etablierten Validierungsprinzipien nennt und dabei Gramatica und den OECD-Leitfaden zitiert, mit dem Zusatz, dass ein kreuzvalidierter Schaetzwert diese Anforderungen nicht erfuellt, wenn die Deskriptorauswahl ausserhalb der Validierungsschleife liegt.");
It("3.2","Geroldinger jetzt auch an der inhaltlich passenden Stelle zitiert.","Die Referenz war in den Conclusions zitiert, aber nicht in dem kurzen Verweisabsatz, der die c-Statistik zusammenfasst. Das ist ergaenzt.");
It("3.3","Literaturliste vollstaendig neu nummeriert und verifiziert.","Jetzt 40 Eintraege, Erstauftritt exakt 1 bis 40, alle zitiert. Die Umnummerierung wurde diesmal mit einer automatischen Zuordnungspruefung durchgefuehrt: fuer mehrere Ankerstellen im Text wurde geprueft, dass die Zitatnummer vor und nach der Umnummerierung auf denselben Literatureintrag zeigt. Bei einem ersten Versuch war Mansouri versehentlich in die Mitte der Liste eingefuegt worden, wodurch alle nachfolgenden Nummern um eins verrutschten; der Fehler wurde durch diese Pruefung entdeckt und die Aenderung von einer gesicherten Fassung aus neu aufgebaut.");

push(H1("4. Nicht geaendert"));
It("4.1","Keine Zahlen, Tabellen oder Abbildungen des Ergebnisteils.","Der GA-Befund steht ausschliesslich im Diskussionsteil und im Supplement. Alle Haupttabellen und alle Abbildungen sind unveraendert.");
It("4.2","Kein vollstaendiges GA-Protokoll gerechnet.","Wie besprochen wurde bewusst nur eine reduzierte Fassung gerechnet, da die vollstaendige Replikation nicht darstellbar ist und das Manuskript sie nicht benoetigt.");

push(H1("5. Offene Punkte"));
It("5.1","Externer Testsatz aus der Literatur.","Offen. Der GA-Befund liefert dafuer ein zusaetzliches Argument: Wo ein externer Satz existiert, ist Selektionsfreiheit unproblematisch, wo keiner existiert, ist sie es nicht.");
It("5.2","Gaucher-Reanalyse als eigenstaendige Mitteilung.","Offen, insbesondere der Richtungsbefund zu m/z 2067.9.");
It("5.3","Weitere Straffung um rund 500 Woerter.","Offen, durch Zusammenlegen der Abschnitte zu Klassenbalance und Klassifikationsmetriken.");
It("5.4","EndNote.","Offen. Jetzt 40 Referenzen in Zitierreihenfolge, feste Nummern statt Feldfunktionen.");

push(new Paragraph({ spacing:{before:260,after:120}, children:[new TextRun({
  text:"Alle Zahlen stammen aus der reproduzierbaren Python-Pipeline. Rohdaten, Skripte und Permutationsausgaben liegen in QSAR_Monoterpene_Gesamtsicherung_2026-08-11.zip.",
  italics:true, size:19, font:FONT})]}));

const doc = new Document({ creator:"Werle et al.", styles:{ default:{ document:{ run:{ font:FONT, size:21 } } } },
  sections:[{ properties:{ page:{ margin:{ top:1440,bottom:1440,left:1440,right:1440 } } }, children:ch }] });
Packer.toBuffer(doc).then(b=>{ fs.writeFileSync("/mnt/user-data/outputs/Werle_JCheminform_Manuskriptaenderungen_2026-08-11.docx", b); console.log("written"); });
