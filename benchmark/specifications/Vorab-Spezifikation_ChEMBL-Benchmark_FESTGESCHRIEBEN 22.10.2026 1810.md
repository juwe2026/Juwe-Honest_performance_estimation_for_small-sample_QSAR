---
title: "Vorab-Spezifikation (vor Kenntnis der Ergebnisse) – ChEMBL-Benchmark"
subtitle: "FESTGESCHRIEBENE FASSUNG vom 22.10.2026, 18:10"
---

# Status und Verfahren

Dies ist die festgeschriebene Fassung. Sie wurde am 22. Oktober 2026 um 18:10 Uhr geschlossen, nachdem die Auswahl der 16 natürlichen Serien abgeschlossen war und bevor das erste Modell gerechnet wurde. Der SHA-256-Wert dieser Datei ist in der Begleitdatei „PRESPEC_SHA256.txt“ festgehalten. Die Fassung wird unverändert als Additional file hinterlegt, mit englischer Übersetzung und dem deutschen Original.

Ab diesem Zeitpunkt wird an den Regeln nichts mehr geändert. Jede spätere Abweichung wird im Abschnitt „Abweichungen“ mit Datum, Grund und Auswirkung protokolliert.

Vor dem Festschreiben war nur Folgendes erlaubt, und nur das ist geschehen: die ChEMBL-Abfrage, das Lesen der Kandidatenartikel zur Ermittlung der Schwelle der Autoren, das Zählen der Klassengrößen und die Berechnung der Gerüst- und Ähnlichkeitskennzahlen für die Auswahl. Kein Modell, kein Filter und keine Permutation wird vorher gerechnet.

Festgelegt von EW, umgesetzt von Claude.

Anlass: Reviewer 1 (Revision 2, Gutachten vom 20.10.2026) verlangt, die allgemeinen Aussagen des Manuskripts an einer größeren Zahl von Datensätzen unterschiedlicher Größe, struktureller Diversität und Aktivitätsart zu prüfen.

# Stand der Entscheidungen

Keine offenen Punkte mehr. Das Dokument kann festgeschrieben werden, sobald die Auswahl der 16 Publikationen abgeschlossen ist.

Entschieden am 21. und 22.10.2026:

- Zugang zu ChEMBL eingerichtet und geprüft (ChEMBL 37, Release 01.05.2026).
- 16 natürliche Serien plus 4 Teilmengen.
- Klassengrenze nach der Schwelle der Autoren, wo vorhanden, sonst Ersatzregel (Ersatzschwellen von Claude festgelegt, von EW übernommen).
- Zahl der Permutationen nach der Regel im Abschnitt „Festgelegte Parameter“ (500, 1.000 oder 2.000).
- Kriterium für H4 wie im Abschnitt „Hypothesen“.
- Die Regel „ein Datensatz je Ziel“ gilt nur für die 16 natürlichen Serien, nicht für die vier Teilmengen.
- Bänder ohne Lücke: Band A 25 bis 35 %, Band B über 35 bis 50 %.
- Je Publikation genau ein Datensatz. Innerhalb der Publikation entscheidet die Hash-Reihenfolge der Assays (siehe „Wahl des Assays innerhalb einer Publikation“).
- Keine Übersichtsarbeiten, keine Arbeiten ohne eigene zusammenhängende Verbindungsserie. Leseliste je Zelle 10 Publikationen.

# Datenquelle

ChEMBL, die zum Abrufzeitpunkt aktuelle Version. Versionsnummer, Abrufdatum und die vollständige Abfrage (SQL bzw. Skript) werden protokolliert und im Code-Archiv hinterlegt.

Einheit eines Datensatzes ist die Kombination aus einem ChEMBL-Dokument (eine Publikation) und einem Assay. Jeder Datensatz stammt also aus einer einzigen Arbeit mit einem einzigen Test. Das entspricht der Situation, die das Manuskript untersucht: eine kleine Serie, ein Labor, ein Test.

# Auswahlregel

## Einschlusskriterien für Kandidaten

- Quelle: Literaturdaten (ChEMBL-Quelle „Scientific Literature“).
- Aktivitätsart (standard_type) aus genau einer der vier Klassen:
  - **MIC:** antibakteriell oder antimykotisch, Assay gegen einen Organismus (target_type ORGANISM).
  - **Enzym:** IC50 gegen ein einzelnes Protein mit Enzymfunktion, confidence_score ≥ 8.
  - **Rezeptor:** Ki gegen ein einzelnes Protein mit Rezeptorfunktion (GPCR, Ionenkanal oder Kernrezeptor), confidence_score ≥ 8.
  - **Zytotoxizität:** IC50, GI50 oder CC50 gegen eine Zelllinie (target_type CELL-LINE).
- Einheitliche Einheit innerhalb des Assays (MIC in µg/mL oder µM, sonst nM).
- Kein data_validity_comment, kein Duplikatvermerk (potential_duplicate = 0).
- Struktur: ChEMBL-Elternstruktur, mit RDKit einlesbar, keine Metallkomplexe, keine Gemische, mindestens 6 Schweratome.
- Mehrfachmessungen derselben Verbindung im selben Assay werden zum Median des Logarithmus zusammengefasst.
- Zensierte Werte (Relation „>“) gelten als nicht-aktiv, wenn die Grenze über der Aktivitätsschwelle liegt, wie bei Monoterpenoiden und Peptiden. Andernfalls wird die Verbindung ausgeschlossen. Werte mit Relation „<“ gelten entsprechend als aktiv, wenn die Grenze unter der Schwelle liegt.
- Größe nach Bereinigung: 30 bis 60 Verbindungen für die natürlichen Serien.
- **Art der Arbeit (beim Lesen zu prüfen):** Zugelassen sind nur Originalarbeiten, deren Aktivitätswerte sämtlich in dieser Arbeit selbst gemessen wurden. Ausgeschlossen sind drei Gruppen:
  1. Übersichtsarbeiten jeder Form, auch als „Review“, „Mini-review“ oder „Perspective“ bezeichnete.
  2. Arbeiten, deren Aktivitätswerte aus anderen Publikationen zusammengetragen sind. Dazu zählen QSAR- und QSSR-Studien, die fremde Datensätze modellieren, und Methoden- oder Assay-Arbeiten, die eine Sammlung bereits publizierter Wirkstoffe vergleichend vermessen.
  3. Arbeiten, deren Verbindungen vor der Messung nach Aktivität vorselektiert wurden, etwa die Treffer aus einer Screening-Bibliothek.
  
  Begründung der ersten beiden Punkte: Das Manuskript untersucht die Auswertung einer Serie, die in einem Labor mit einem Assay gemessen wurde. Eine zusammengetragene Werteliste erfüllt das nicht. Ein gemeinsames Grundgerüst wird dagegen nicht verlangt, denn chemische Vielfalt ist ein Faktor des Versuchsplans, und auch der Gaucher-Datensatz des Manuskripts hat keine gemeinsame Chemie.
  
  Begründung des dritten Punktes: Bei Screening-Treffern liegt vor unserer Klassendefinition bereits eine labelabhängige Auswahl der Verbindungen, die wir weder kontrollieren noch beschreiben können. In einer Arbeit über labelabhängige Auswahl wäre das besonders unglücklich.

## Klassenverteilung

- Die kleinere Klasse (aktiv oder nicht-aktiv) umfasst mindestens 25 % und höchstens 50 % der Verbindungen und mindestens 10 Verbindungen.
- Die beiden Bänder schließen lückenlos aneinander an: Band A 25 bis 35 %, Band B über 35 bis 50 %. Jeder zulässige Assay fällt damit in genau ein Band.
- Begründung: Im Manuskript war P. aeruginosa mit 3 von 31 aktiven (10 %) nicht auswertbar und lief nur mit. Die knappsten auswertbaren Endpunkte hatten 5 und 6 aktive. Die Grenze liegt bewusst deutlich darüber.
- Die Schwelle wird nicht verschoben, um einen Datensatz passend zu machen. Ein Datensatz, der die Grenze verfehlt, wird mit diesem Grund ausgeschlossen.

## Strukturelle Diversität

Berechnet werden die mittlere paarweise Tanimoto-Ähnlichkeit (Morgan, Radius 2, 2.048 Bit) und der Anteil der Verbindungen im häufigsten generischen Bemis-Murcko-Gerüst.

- **Kongener:** häufigstes generisches Gerüst ≥ 50 % der Verbindungen.
- **Divers:** häufigstes generisches Gerüst < 25 % und mittlere Tanimoto-Ähnlichkeit < 0,30.
- Datensätze dazwischen werden nicht verwendet.

## Versuchsplan

Ziel sind genau 16 natürliche Serien und 4 Teilmengen, zusammen 20 Datensätze. Vier Aktivitätsarten × zwei Diversitätsstufen × zwei Bänder der Klassenverteilung ergeben 16 Zellen mit je einem Datensatz:

| Faktor | Stufen |
|---|---|
| Aktivitätsart | MIC, Enzym-IC50, Rezeptor-Ki, Zytotoxizität |
| Diversität | kongener, divers |
| Anteil der kleineren Klasse | Band A 25 bis 35 %, Band B über 35 bis 50 % |

Zusätzliche Einschränkungen: höchstens ein Datensatz je Dokument, und höchstens ein Datensatz je Zielprotein bzw. je Organismus bzw. je Zelllinie.

Die Zielregel wird in zwei Durchgängen angewandt. Im **ersten Durchgang** gilt sie strikt. Bleiben danach Zellen leer, folgt ein **zweiter Durchgang**, in dem sich ein Ziel wiederholen darf, aber höchstens einmal, also höchstens zwei Datensätze je Ziel insgesamt. Der zweite Durchgang geht die Publikationen wieder in der Hash-Reihenfolge durch und berücksichtigt nur noch leere Zellen.

Begründung: Die Zielregel dient der biologischen Breite des Versuchsplans. Die statistische Unabhängigkeit der Datensätze sichert dagegen die Regel „ein Datensatz je Publikation“, denn zwei Arbeiten zum selben Ziel haben verschiedene Verbindungen, verschiedene Deskriptormatrizen und stammen aus verschiedenen Laboren. Ein vollständig besetzter Versuchsplan ist für die Aussage wichtiger als die strikte Einmaligkeit des Ziels. Die Auffüllregel des nächsten Abschnitts greift erst, wenn auch der zweite Durchgang eine Zelle nicht besetzen kann. Beide Einschränkungen gelten nur für die 16 natürlichen Serien. Die vier Teilmengen der Größenvariation sind davon ausgenommen, weil sie explorativ ausgewertet werden und nicht in die Tests von H1 bis H4 eingehen. Ein Elterndatensatz darf also dasselbe Ziel oder denselben Organismus betreffen wie eine der natürlichen Serien, aber nicht dasselbe Dokument.

## Wahl des Assays innerhalb einer Publikation

Eine Publikation enthält fast immer mehrere Assays, bei MIC-Serien meist vier bis sechs Keime, bei Zytotoxizität mehrere Zelllinien. Alle Assays einer Publikation beruhen auf **denselben Verbindungen**. Sie haben deshalb dieselbe Deskriptormatrix, dieselbe Korrelationsgruppierung und dasselbe Redundanzverhältnis G/p. Nur der Labelvektor unterscheidet sie.

Daraus folgt: **Aus einer Publikation geht genau ein Datensatz in die Auswertung ein.** Zwei Assays derselben Arbeit wären für die Prüfung von H3 keine unabhängigen Beobachtungen, und die Tests über die Datensätze hinweg setzen unabhängige Blöcke voraus.

Welcher Assay das ist, entscheidet folgende Regel. Die Assays einer Publikation werden in der Hash-Reihenfolge durchgegangen. Genommen wird der erste Assay, der alle Auswahlkriterien erfüllt, einschließlich der Klassenverteilung, **und** in ein Band fällt, das in der betreffenden Zelle noch frei ist. Gibt es keinen solchen Assay, wird die Publikation nicht verwendet.

Begründung: Der Versuchsplan ist nur aussagekräftig, wenn seine Zellen besetzt sind. In den diversen Zellen stehen je nach Aktivitätsart nur 27 bis 82 Kandidaten zur Verfügung, und ein erheblicher Teil der Assays verfehlt die Klassengrenze. Eine Regel, die eine Publikation schon dann verwirft, wenn ihr erster geeigneter Assay in ein bereits besetztes Band fällt, würde leere Zellen in Kauf nehmen. Ein unvollständiger Plan ist für die Aussage schädlicher als die Auswahl des Assays nach einem Merkmal, das ohnehin ein Faktor des Plans ist.

Die Auswahl verwendet ausschließlich Merkmale, die vor jeder Modellrechnung feststehen: Aktivitätsart, Struktur, Größe, Diversität und Klassenverteilung. Kein Ergebnis einer Modellrechnung, also weder Q² noch AUC noch Trefferzahl oder Inflationsfaktor, geht in die Auswahl ein. Diese Größen sind zum Zeitpunkt der Auswahl nicht berechnet.

Im Blatt „Selection“ wird **jeder geprüfte Assay** aufgeführt, mit Klassenanteil, Band und Grund des Übergehens, sowie je gewähltem Datensatz die Angabe, der wievielte geeignete Assay seiner Publikation er war.

## Reihenfolge der Publikationen (gegen Rosinenpickerei)

Zwei der Kriterien stehen erst beim Lesen des Artikels fest: die Schwelle der Autoren und die Art der Arbeit. Die Reihenfolge entscheidet deshalb, in welcher Folge gelesen wird.

Alle Kandidaten werden in eine feste Pseudozufallsreihenfolge gebracht: SHA-256 von „QSAR-R2-“ plus ChEMBL-Dokument-ID plus Assay-ID, aufsteigend sortiert. Der Rang einer Publikation ist der Rang ihres ersten Assays. In dieser Reihenfolge werden die Publikationen geprüft, jede nach der Regel des vorigen Abschnitts. Die Prüfung endet, sobald alle 16 Plätze besetzt sind. Weitere Publikationen werden nicht mehr geprüft.

Auffüllregel: Bleibt eine Zelle auch nach dem zweiten Durchgang leer, erhält sie den nächsten passenden Kandidaten derselben Aktivitätsart in der Hash-Reihenfolge, ohne Rücksicht auf Diversität und Band. Gibt es auch den nicht, wird der nächste passende Kandidat einer beliebigen Aktivitätsart genommen. Jede Auffüllung wird mit Grund protokolliert. Am Ende stehen immer 16 natürliche Serien.

„Passend“ wird ausschließlich anhand der Auswahlkriterien entschieden, vor jeder Modellrechnung. Ein Datensatz wird nach Beginn der Rechnungen nur aus technischen Gründen ersetzt (zum Beispiel Deskriptorberechnung für mehr als 10 % der Verbindungen fehlgeschlagen), dann durch den nächsten passenden Kandidaten derselben Zelle. Ein Ersatz wegen eines Ergebnisses ist ausgeschlossen.

Alle Kandidaten mit Ein- oder Ausschluss und Grund werden im Blatt „Selection“ der Ergebnisdatei aufgeführt.

## Größenvariation (4 Teilmengen, explorativ)

- Zwei Elterndatensätze: der größte zulässige MIC-Datensatz und der größte zulässige Enzym-IC50-Datensatz mit je mindestens 100 Verbindungen und mindestens 30 % in der kleineren Klasse. Sie dürfen nicht aus einem Dokument stammen, das schon eine der 16 natürlichen Serien stellt. Dasselbe Ziel oder derselben Organismus ist dagegen zulässig (siehe „Versuchsplan“).
- Aus jedem Elterndatensatz werden zwei Teilmengen gezogen, n = 25 und n = 60, geschichtet nach Klasse mit 30 % in der kleineren Klasse (8 bzw. 18 Verbindungen). Bei n = 25 wird die Untergrenze von 10 Verbindungen bewusst unterschritten. Das ist der einzige Fall.
- Ziehung mit festem Seed (siehe „Seeds“). Die Teilmengen dienen nur der explorativen Auswertung des Größeneffekts und gehen nicht in die Tests von H1 bis H4 ein.

# Klassendefinition

Die Klassengrenze wird je Datensatz in einer festen Rangfolge bestimmt, vor jeder Modellrechnung.

**Stufe 1 – Schwelle der Autoren.** Nennt der Artikel ausdrücklich eine Grenze für aktiv und nicht-aktiv, wird sie übernommen. Das gilt auch für eine Schwelle, die die Autoren erst in der Studie ermittelt haben, und für eine erkennbare Einteilung in der Darstellung (zum Beispiel „active“ gegenüber „inactive“ in einer Tabelle). Im Blatt „Selection“ werden die Textstelle wörtlich, die Seite bzw. Tabelle und die daraus abgeleitete Regel festgehalten. Ob die Autoren die Grenze vorab oder aus den Daten festgelegt haben, ist für diese Studie unerheblich. Entscheidend ist, dass sie unabhängig von unseren Modellen feststeht.

**Stufe 2 – Ersatzregel.** Nennt der Artikel keine Grenze, gilt die übliche Grenze der Aktivitätsart:

| Aktivitätsart | aktiv | Begründung |
|---|---|---|
| MIC | ≤ 8 µg/mL (µM-Werte über die Molmasse umgerechnet) | übliche Screening-Grenze, wie bei den Peptiden (Additional file 5) |
| Enzym-IC50 | pChEMBL ≥ 6 (≤ 1 µM) | übliche Grenze für Wirksamkeit in der Leitstrukturphase |
| Rezeptor-Ki | pChEMBL ≥ 7 (≤ 100 nM) | Rezeptor-Ki-Werte liegen typischerweise eine Größenordnung niedriger |
| Zytotoxizität | pChEMBL ≥ 5 (≤ 10 µM) | übliche Grenze in Zellassays |

Diese Schwellen sind Screening-Grenzen ohne klinische Bedeutung, wie im Manuskript. Dasselbe Vorgehen wurde bereits im Manuskript angewandt: bei den Monoterpenoiden die Grenze aus Ref. 14, bei den Peptiden mangels Grenze der Autoren eine vorab festgelegte (Additional file 5).

Die Herkunft der Grenze (Stufe 1 oder 2) wird je Datensatz als Merkmal geführt.

In beiden Stufen gilt die Untergrenze der Klassenverteilung. Ergibt die Grenze weniger als 25 % oder weniger als 10 Verbindungen in der kleineren Klasse, wird der Kandidat ausgeschlossen. Die Grenze wird nicht verschoben.

**Literaturbeschaffung:** Artikel mit freiem Volltext liest Claude selbst. Die übrigen werden von den Autoren beschafft und in den Ordner „Revision 2“ gelegt. Ein Kandidat wird nicht wegen fehlenden Zugangs übersprungen.

**Umfang der Leseliste:** Je Zelle werden die ersten 10 Publikationen der Hash-Reihenfolge vorgelegt. Grund: Die ersten gelesenen Artikel zeigen, dass ein erheblicher Teil der Kandidaten erst beim Lesen ausscheidet, vor allem als Übersichtsarbeit oder Methodenarbeit. Gelesen wird in der Rangfolge, bis die Zelle besetzt ist. Die übrigen Artikel der Zelle müssen dann nicht mehr gelesen werden.

**Teilmengen der Größenvariation:** Für die beiden Elterndatensätze gilt dieselbe Rangfolge.

**Sensitivitätsanalyse:** Median-Regel (aktiv = besser als der Median des Datensatzes). Sie ergibt etwa 50 % aktive und wird für alle 16 Datensätze gerechnet, aber nur deskriptiv berichtet.

# Deskriptoren

- Wie beim zweiten MIC-Datensatz: Mordred (nur 2D) plus RDKit, aus der ChEMBL-Elternstruktur als kanonischem SMILES. Keine Protonierungs- oder pH-Korrektur, keine 3D-Deskriptoren.
- Bereinigung wie im Manuskript: nicht-numerische Spalten entfernen, Spalten mit mehr als 50 % fehlenden Werten entfernen, Rest median-imputieren, konstante Spalten entfernen.
- Zusätzlich Morgan-Fingerprints (Radius 2, 2.048 Bit) für die Random-Forest-Kontrolle.
- Je Datensatz protokolliert: Zahl der Deskriptoren roh und bereinigt, n:p, Zahl der Korrelationsgruppen, Redundanzverhältnis G/p (Gruppen geteilt durch bereinigte Deskriptoren).

# Zu rechnende Analysen je Datensatz

Identisch mit Additional file 5, ohne Route B (keine Literaturdeskriptoren):

- E1 Route A1 (Gruppierung zuerst, Repräsentant nach Effektgröße, dann univariater Filter), scheinbar und genestet
- E1 Route A1u (Repräsentant nach Zentralität, labelfrei), scheinbar und genestet
- E1 Route A2 (Filter zuerst, dann Gruppierung), scheinbar und genestet
- E2 Ohne jede Auswahl: PLS-DA auf allen bereinigten Deskriptoren
- E2 Ohne Filter: PLS-DA auf den labelfreien Gruppenrepräsentanten
- E4 Random Forest (300 Bäume, min_samples_leaf = 1) auf Morgan-Fingerprints, ohne Auswahl, OOB
- E6 Nullverteilung der univariaten Trefferzahl für A2, A1 und A1u, Inflationsfaktoren nominal (5 %) und auf die eigene Null-Trefferrate kalibriert

Kennzahlen: Q² und AUC (für Tests), dazu Sensitivität, Spezifität, balancierte Genauigkeit und MCC bei Schwelle 0,5, wie im Manuskript.

# Festgelegte Parameter

- PLS-DA mit 1 latenter Variable, Leave-one-out, Autoskalierung nur auf den Trainingsverbindungen jedes Folds, Spalten mit Trainings-SD < 1e-6 im Fold verworfen (labelfrei).
- Univariater Filter: Mann-Whitney p ≤ 0,05 und |r_rb| > 0,3, in jedem Fold neu.
- Korrelationsgruppierung: |Pearson| oder |Spearman| > 0,7, most-connected-first, Gleichstand → größte Varianz.
- Folds ohne Deskriptor oder mit nur einer Klasse: Vorhersage 0,5 (nicht der Trainingsmittelwert), wie im Manuskript.
- Permutationen: scheinbare Modelle 1.000, Random Forest 500, Trefferzahl 1.000. Genestete Modelle nach folgender Regel:
  - Vor dem ersten Rechenlauf wird die Laufzeit gemessen, und zwar ausschließlich mit permutierten Labels am Datensatz mit dem größten Deskriptorpool. Dabei wird kein Ergebnis mit den beobachteten Labels berechnet.
  - Daraus wird die Gesamtrechenzeit für alle 20 Datensätze hochgerechnet.
  - Die Hochrechnung wird vor der Entscheidung berichtet, getrennt für 500, 1.000 und 2.000 Permutationen. EW entscheidet daraufhin, welche Zahl verwendet wird. Richtwert ist die größte Zahl, deren hochgerechnete Gesamtzeit 48 Stunden nicht übersteigt. Stehen mehr Zeit und Rechenleistung zur Verfügung, kann auch eine höhere Zahl gewählt werden. Die Zahl gilt einheitlich für alle Datensätze.
  - Die Entscheidung stützt sich ausschließlich auf Laufzeiten, die mit permutierten Labels gemessen wurden. Kein Ergebnis mit den beobachteten Labels liegt zu diesem Zeitpunkt vor.
  - Monte-Carlo-Genauigkeit eines p-Werts um 0,05 (halbe Breite des 95-%-Intervalls): 0,019 bei 500, 0,014 bei 1.000, 0,010 bei 2.000 Permutationen. Kleinster erreichbarer p-Wert: 0,002, 0,001 bzw. 0,0005.
  - Messwerte, Hochrechnung und gewählte Zahl werden im Abschnitt „Abweichungen“ bzw. im README protokolliert.
- Paarvergleiche der Q²: gepaarter, nach Klasse geschichteter Bootstrap mit 20.000 Ziehungen, 95-%-Perzentilintervall, wie im Manuskript.

# Seeds

numpy default_rng([93, Datensatznummer, Analysecode, Permutationsindex]).
Datensatznummern 1 bis 16 (natürliche Serien, in der Reihenfolge des Versuchsplans), 21 bis 24 (Teilmengen).
Analysecodes: A1 = 1, A1u = 2, A2 = 3, ohne Auswahl = 4, ohne Filter = 5, RF = 6, Trefferzahl = 7, Bootstrap = 8, Teilmengenziehung = 9.
Sensitivitätsanalyse Median-Regel: Analysecode plus 100.

# Hypothesen

Primäre Analysemenge sind die 16 natürlichen Serien. Die drei Datensätze des Manuskripts (Monoterpenoide, Peptide, Gaucher) gehen nicht in die Tests ein. Sie werden in den Abbildungen markiert. Die Hypothesen sind die zentralen Aussagen des Manuskripts. Ihre Prüfung ist eine externe Replikation. Alle Tests sind einseitig in der vorhergesagten Richtung, α = 0,05, Holm-Korrektur über die vier primären Tests H1 bis H4.

**H1 – Optimismus durch Selektion.** Unter Route A1 liegt die scheinbare Q² über der genesteten.
Test: Wilcoxon-Vorzeichen-Rang-Test über die 16 Differenzen (scheinbar minus genestet). Sekundär dasselbe für A1u und A2 sowie die Zahl der Datensätze mit positiver Differenz.

**H2 – Reihenfolge der Inflation.** Der kalibrierte Inflationsfaktor der Trefferzahl folgt der Ordnung A2 > A1 > A1u.
Test: Page-Trendtest für geordnete Alternativen mit den 16 Datensätzen als Blöcken. Sekundär: paarweise Wilcoxon-Tests A2 gegen A1 und A1 gegen A1u, Zahl der Datensätze mit vollständig vorhergesagter Ordnung.

**H3 – Redundanz bestimmt die Stärke.** Je redundanter der Deskriptorpool (kleineres G/p), desto größer der kalibrierte Inflationsfaktor unter Route A2.
Test: Spearman-Korrelation zwischen G/p und dem Logarithmus des Inflationsfaktors über die 16 Datensätze, vorhergesagt negativ. Sekundär dasselbe für A1 und A1u. Die drei Datensätze des Manuskripts werden in der Abbildung als Vorhersage eingetragen, die sich bewähren muss.

**H4 – Kein Gewinn durch überwachte Selektion.** Keine gefilterte Route übertrifft die PLS-DA auf den labelfreien Gruppenrepräsentanten (Kontrolle ohne Filter) über die Stichprobenunsicherheit hinaus.
Zählung: Datensätze, in denen für mindestens eine der drei Routen das 95-%-Bootstrapintervall der Differenz der genesteten Q² (Route minus Kontrolle) vollständig über 0 liegt.
Entscheidungskriterium: H4 gilt als widerlegt, wenn dies in mindestens 4 der 16 Datensätze eintritt. Begründung: Ein zweiseitiges 95-%-Intervall liegt unter H4 mit höchstens 2,5 % Wahrscheinlichkeit vollständig über 0. Bei drei Routen je Datensatz sind es im ungünstigsten Fall (unabhängige Routen) 7,5 % je Datensatz, also im Mittel 1,2 von 16 Datensätzen. Die Wahrscheinlichkeit für 4 oder mehr beträgt dann 0,028 (Binomialverteilung), liegt also unter 0,05.
Sekundär: dieselbe Zählung gegen die PLS-DA ohne jede Auswahl und gegen den Random Forest. Zusätzlich als Effektschätzung der Median der Differenzen (beste gefilterte Route minus Kontrolle) über die 16 Datensätze mit Bootstrap-Intervall über die Datensätze.
Lesart: Wird H4 nicht widerlegt, ist das kein Beweis der Gleichwertigkeit. Bei rund 40 Verbindungen sind die Intervalle breit. Aussagekräftig ist deshalb vor allem die Effektschätzung: liegt sie nahe 0 oder darunter, bringt die Selektion im Mittel nichts.

**Explorativ (ohne Test-Anspruch):**

- Einfluss des Bandes der Klassenverteilung (Band A gegen Band B) auf Optimismus und Inflation.
- Einfluss der Diversität und der Aktivitätsart.
- Größeneffekt aus den vier Teilmengen (n = 25 gegen 60 innerhalb desselben Elterndatensatzes).
- Ob unter A1u ein Überschuss an Treffern gegenüber der labelfreien Nullverteilung mit einer nominal signifikanten genesteten Q² zusammenfällt, wie bei S. aureus und S. pneumoniae.
- Herkunft der Grenze (Autoren gegenüber Ersatzregel) als erklärendes Merkmal.
- Median-Regel als Sensitivitätsanalyse.

# Berichtspflichten

- Alle Datensätze und alle Analysen werden berichtet, auch wenn eine Hypothese scheitert.
- Abweichungen von dieser Spezifikation werden mit Datum, Grund und Auswirkung im Abschnitt „Abweichungen“ der Ergebnisdatei protokolliert.
- Ergebnisdatei: eine Excel-Arbeitsmappe mit den Blättern README, Selection, Summary, All_results_long und DS01 bis DS16 sowie DS21 bis DS24 im identischen Layout. Das Blatt „Selection“ enthält das vollständige Auswahlprotokoll einschließlich aller übergangenen Assays. SMILES, ChEMBL-IDs, Rohwerte, Klassen, Deskriptormatrizen und Skripte gehen ins Code-Archiv bzw. nach Zenodo.

# Abweichungen

(leer – wird nach dem Festschreiben geführt)
