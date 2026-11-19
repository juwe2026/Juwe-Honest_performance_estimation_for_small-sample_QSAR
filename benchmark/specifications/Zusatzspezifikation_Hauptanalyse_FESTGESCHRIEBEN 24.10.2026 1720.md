---
title: "Zusatzspezifikation – Hauptanalyse ohne Kontroll- und Referenzsubstanzen"
subtitle: "FESTGESCHRIEBENE FASSUNG vom 24.10.2026, 17:20"
---

# Status

Diese Fassung wurde am 24.10.2026 um 17:20 Uhr festgeschrieben, bevor ein Modell dieser Analyse gerechnet wurde. Ihr SHA-256-Wert steht in „ZUSATZSPEZ_SHA256.txt". Festgelegt von EW, umgesetzt von Claude.

Die Vorab-Spezifikation vom 22.10.2026, 18:10 (SHA-256 a13bb2c0b513289b160de744bd8d6bc47662d46553ab5a2fe3db6d336a0ae419) bleibt unverändert hinterlegt. Diese Zusatzspezifikation ergänzt sie um eine Regel, die dort fehlte. Alle übrigen Festlegungen der Vorab-Spezifikation gelten unverändert: Auswahlregel, Klassendefinition, Parameter, Hypothesen H1 bis H4 und ihre Tests.

**Entscheidung EW (24.10.2026, 14:23): Die Analyse nach dieser Zusatzspezifikation ist die Hauptanalyse des Manuskripts.** Die Analyse nach der Vorab-Spezifikation allein wird vollständig im Supplement berichtet. Im Methodenteil steht ein Satz zu dieser Ergänzung mit Verweis auf das Protokoll. Die Gutachter werden kurz informiert. Diese Entscheidung fiel vor jeder Rechnung dieser Analyse und gilt unabhängig von ihrem Ergebnis.

**Anlass und Zeitpunkt.** Die Lücke wurde nach Kenntnis der Ergebnisse der Vorab-Analyse erkannt. Anlass war die Untersuchung numerisch entgleister Vorhersagen. Dabei zeigte sich, dass mehrere Serien Arzneistoffe enthalten, die die Autoren nur zum Vergleich oder zur Absicherung des Tests mitgemessen haben. Die Ergänzung betrifft nur, welche Verbindungen zu einer Serie gehören. Kein Ergebnis einer Modellrechnung geht in die Regel oder in die Auswahl ein.

# Regel: Was gehört zur untersuchten Serie?

**Zur Serie gehören**

1. alle Verbindungen, die die Arbeit als neue Verbindungen der Serie vorstellt,
2. die Ausgangs-, Leit- oder Vorläuferstrukturen, von denen die Serie nach Angabe des Artikels abgeleitet ist, und bekannte Vertreter derselben Substanzklasse, die in der Arbeit im selben Test zusammen mit den neuen Verbindungen gemessen wurden (Ergänzung EW, 24.10.2026, 08:45).

**Nicht zur Serie gehören**

3. Arzneistoffe, klinische Kandidaten oder Werkzeugsubstanzen anderer Substanzklassen, die zum Vergleich, als Referenz, als Standard oder als Positivkontrolle mitgemessen wurden, zum Beispiel zur Charakterisierung von Resistenzmustern oder zur Absicherung des Tests, gleich ob gekauft oder nachsynthetisiert,
4. Verbindungen, deren Messwert der Artikel aus einer anderen Publikation übernimmt.

**Entscheidend ist die Herkunft und die Zugehörigkeit zur Substanzklasse, nicht das Wort.** Entschieden wird am Artikeltext, mit Textstelle im Protokoll. Als Hilfe dienen die ChEMBL-Vorzugsnamen, die Verbindungsschlüssel des Artikels und die strukturelle Ähnlichkeit zur übrigen Serie (Morgan, Radius 2, 2.048 Bit, Tanimoto).

# Anwendung

## Entscheidungen je Publikation (geprüft an den Verlags-PDFs, von EW am 24.10.2026, 14:13 und 14:23 bestätigt)

| Publikation | ausgeschlossen | Punkt |
|---|---|---|
| Noji 2020 (JAK3) | keine | |
| Ragab 2016 (COX-2) | Celecoxib, Diclofenac | 3 |
| Garofalo 2012 (VEGFR2) | keine | |
| Rewcastle 1995 (EGFR) | keine | |
| Brown 2011 (M. tuberculosis) | Isoniazid, Ethambutol | 3 |
| Hansen 2019 (E. coli) | 12 Antibiotika (CIP, MEM, PIP, IMP, ATM, CAZ, AMP, AK, FEP, TGC, TE, CTX) | 3 |
| Kozic 2012 (M. tuberculosis) | Isoniazid | 3 |
| Zhao 2023 (C. albicans) | Fluconazol | 3 |
| Bonifazi 2023 (D3) | Oliceridin (3), TRV734 (5) | 3 |
| Bonifazi 2021 (MOR) | Fentanyl, Haloperidol (5). Loperamid (3), Diphenoxylat (4) und Methadon (6) bleiben nach Punkt 2 | 3 |
| Kai 2007 (CB2) | WIN55212-2 | 3 |
| Haadsma-Svensson 1995 (D2) | 8-OH-DPAT | 3 |
| Jardim 2015 (HCT-116) | Doxorubicin | 3 |
| Díaz-Gavilán 2008 (MCF7) | 5-Fluorouracil | 3 und 4 |
| Lu 2010 (RH7777) | keine | |
| Kamal 2015 (HeLa) | Doxorubicin | 3 |
| Ogata 1987 (A. fumigatus), nachgerückt | Clotrimazol, Miconazol, Econazol, Croconazol („known imidazole antimycotics", Tabelle I) | 3 |
| Suckling 2007 (5-HT7), nachgerückt | keine | |
| Spink 2015 (S. aureus), Elterndatensatz MIC | Vancomycin, Linezolid (Tabelle 1, Positivkontrolle im Tierversuch) | 3 |
| Deau 2023 (CLK3), Elterndatensatz Enzym | Leucettine L90, L91, L92 („Data for Leucettines L90−L92 are from ref 65"). Die Referenz-Kinasehemmer der Selektivitätstabelle sind im CLK3-Assay nicht enthalten | 4 |
| Carradori 2017, bisheriger Elterndatensatz MIC | Fluconazol, Ketoconazol, Clotrimazol, Miconazol, Tioconazol, Amphotericin B, Ceftazidim | 3 |

Grenzfälle, die nach Punkt 1 oder 2 in der Serie bleiben: Delgocitinib (Noji 2020), Isoxyl (43, Brown 2011), Normethadon (19, Bonifazi 2021), Leitverbindung 1 (U93385, Haadsma-Svensson 1995), β-, α- und nor-β-Lapachon (Jardim 2015), 1a und 1b (Lu 2010).

## Mechanische Wiederholung der Auswahl

Die Auswahl der Vorab-Spezifikation wird mit dieser Regel vollständig wiederholt: dieselbe Hash-Reihenfolge, dieselben zwei Durchgänge, dieselbe Auffüllregel, dieselben Größen-, Klassen- und Diversitätskriterien. Für jeden geprüften Assay werden zuerst die Fremdsubstanzen entfernt, dann Größe, Klassen und Diversität neu berechnet. Für noch nicht gelesene Publikationen dient als Suchkriterium: Verbindung mit ChEMBL-Vorzugsname und höchster Ähnlichkeit zu den unbenannten Verbindungen des Assays unter 0,60. Wird ein so bereinigter Assay gewählt oder ändert er dadurch seinen Status, wird der Artikel gelesen. Skript: auswahl_neu.py, Protokoll: kontrollen/auswahl_neu.json.

**Ergebnis:**

- 14 der 16 Zellen behalten ihre Publikation und ihren Assay.
- Zelle MIC / divers / Band B: Hansen 2019 hat ohne die Antibiotika 18 Verbindungen und ist unzulässig. Nachgerückt ist Ogata 1987 (Assay CHEMBL649170, A. fumigatus). Mit den vier Azolen hatte die Arbeit 62 Verbindungen und war zu groß. Ohne sie: 58 Verbindungen, 22 aktiv, 36 nicht aktiv, divers.
- Zelle Rezeptor / divers / Band A: Bonifazi 2023 hat ohne Oliceridin und TRV734 eine mittlere Tanimoto-Ähnlichkeit von 0,308 und fällt aus „divers" heraus. Nachgerückt ist Suckling 2007 (Assay CHEMBL886703, 5-HT7), der nächste Kandidat der Hash-Reihenfolge: 32 Verbindungen, 10 aktiv, 22 nicht aktiv.
- Elterndatensatz MIC: Carradori 2017 hat ohne Referenzsubstanzen 116 Verbindungen. Größter zulässiger MIC-Datensatz ist jetzt Spink 2015 (Assay CHEMBL3414284) mit 119 Verbindungen und 34,5 % in der kleineren Klasse.
- Elterndatensatz Enzym: unverändert Deau 2023 (Assay CHEMBL5363956, CLK3), ohne L90 bis L92 mit 229 Verbindungen.

Klassendefinition: In allen Publikationen außer Spink 2015 gilt die Ersatzregel. Spink 2015 nennt „Active compounds (MIC ≤ 8 μg/mL)". Das ist Stufe 1 und fällt mit der Ersatzregel zusammen.

## Datensätze der Hauptanalyse (datasets_k/)

| Datensatz | Publikation | n | aktiv / nicht aktiv | Fremdsubstanzen entfernt |
|---|---|---|---|---|
| DS01 | Noji 2020 | 50 | 34 / 16 | 0 (übernommen) |
| DS02 | Ragab 2016 | 34 | 14 / 20 | 2 |
| DS03 | Garofalo 2012 | 45 | 12 / 33 | 0 (übernommen) |
| DS04 | Rewcastle 1995 | 40 | 23 / 17 | 0 (übernommen) |
| DS05 | Brown 2011 | 42 | 30 / 12 | 2 |
| DS06 | Ogata 1987 (neu) | 58 | 22 / 36 | 4 |
| DS07 | Kozic 2012 | 42 | 13 / 29 | 1 |
| DS08 | Zhao 2023 | 36 | 15 / 21 | 1 |
| DS09 | Suckling 2007 (neu) | 32 | 10 / 22 | 0 |
| DS10 | Bonifazi 2021 | 45 | 19 / 26 | 2 |
| DS11 | Kai 2007 | 34 | 10 / 24 | 1 |
| DS12 | Haadsma-Svensson 1995 | 30 | 11 / 19 | 1 |
| DS13 | Jardim 2015 | 34 | 10 / 24 | 1 |
| DS14 | Díaz-Gavilán 2008 | 48 | 22 / 26 | 1 |
| DS15 | Lu 2010 | 39 | 27 / 12 | 0 (übernommen) |
| DS16 | Kamal 2015 | 32 | 13 / 19 | 1 |
| DS21, DS22 | Spink 2015, Teilmengen n = 25 und 60 | 25, 60 | 8 / 17, 18 / 42 | Eltern: 2 |
| DS23, DS24 | Deau 2023, Teilmengen n = 25 und 60 | 25, 60 | 17 / 8, 42 / 18 | Eltern: 3 |

# Rechnung

- Pipeline, Parameter, Permutationszahlen (genestete Modelle 2.000, scheinbare Modelle und Trefferzahl 1.000, Random Forest 500) und Bootstrap (20.000) wie in der Vorab-Spezifikation. Bereinigung, Korrelationsgruppierung und Repräsentanten werden für jeden veränderten Datensatz neu bestimmt.
- Seeds: numpy default_rng([93, Datensatznummer, Analysecode + 200, Permutationsindex]). Bootstrap: Analysecode 8 + 200 = 208. Teilmengenziehung: 9 + 200 = 209.
- DS01, DS03, DS04 und DS15 sind identisch mit der Vorab-Analyse. Ihre Ergebnisse werden übernommen, nicht neu gerechnet.
- Neu gerechnet werden DS02, DS05 bis DS14, DS16 und DS21 bis DS24.

# Auswertung

- H1 bis H4 mit denselben Tests, Kennzahlen und der Holm-Korrektur wie in der Vorab-Spezifikation, auf den 16 Serien. H4-Grenze 4 von 16.
- Robustheit wie für die Vorab-Analyse beschlossen (EW, 24.10.2026, 04:59): AUC statt Q² für H1 und H4 sowie Ausschluss der Serien mit entgleisten Vorhersagen (Kriterien R1 und R2).
- Die Vorab-Analyse wird im Supplement vollständig berichtet, mit Gegenüberstellung je Hypothese.

# Laufzeit

Hochrechnung aus den gemessenen Paketzeiten der Vorab-Rechnung, skaliert mit n² und der Deskriptorzahl: etwa 23 Stunden, mit Zuschlag 26 Stunden, also etwa 14 Blöcke zu 110 Minuten. Start nur nach Freigabe durch EW.

# Abweichungen

(leer – wird nach dem Festschreiben geführt)
