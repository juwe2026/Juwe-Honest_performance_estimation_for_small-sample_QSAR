---
title: "Ergänzung zur Spezifikation – Sensitivitätsanalyse Median-Regel"
subtitle: "FESTGESCHRIEBENE FASSUNG vom 25.10.2026, 21:39"
---

# Status

Diese Fassung wurde am 25.10.2026 um 21:39 Uhr festgeschrieben, nach Freigabe durch EW (25.10.2026) und bevor ein Modell dieser Analyse gerechnet wurde. Ihr SHA-256-Wert steht in „ERGAENZUNG_MEDIAN_SHA256 25.10.2026 2139.txt“. Der Entwurf vom 25.10.2026 21:26 wurde unverändert übernommen, nur dieser Status wurde angepasst. Die Klassen nach der Median-Regel sind bestimmt, weil dafür nur die Messwerte nötig sind. Gerechnet wurden außerdem nur Laufzeiten mit permutierten Klassen. Kein Modell mit den Klassen der Median-Regel wurde gerechnet.

Grundlage ist die Vorab-Spezifikation vom 22.10.2026 (SHA-256 a13bb2c0…ae419), Abschnitt „Klassendefinition": „Sensitivitätsanalyse: Median-Regel (aktiv = besser als der Median des Datensatzes). Sie ergibt etwa 50 % aktive und wird für alle 16 Datensätze gerechnet, aber nur deskriptiv berichtet." Seeds: „Sensitivitätsanalyse Median-Regel: Analysecode plus 100." Die Vorab-Spezifikation legt nicht fest, wie zensierte Werte und Werte genau auf dem Median zu behandeln sind. Das regelt diese Ergänzung.

Angewandt wird die Median-Regel auf die 16 Serien der Hauptanalyse (Zusatzspezifikation vom 24.10.2026, SHA-256 36668775…df87).

Festgelegt von EW, umgesetzt von Claude.

# Regel

1. **Datensätze.** Dieselben Verbindungen und dieselben Deskriptormatrizen wie in der Hauptanalyse. Nur die Klassen ändern sich. So zeigt der Vergleich allein die Wirkung der Klassengrenze.
2. **Grenze.** Median der erfassten Aktivitätswerte der Serie. Je Verbindung gilt der Wert wie in der Hauptanalyse: Median der logarithmierten Einzelwerte, bei ausschließlich zensierten Angaben die Grenze. Zensierte Werte gehen mit ihrer Grenze in den Median ein.
3. **Aktiv** ist eine Verbindung, deren Wert echt kleiner als der Median ist (niedriger Wert = stärkere Wirkung bei MIC, IC50 und Ki).
4. **Gleich dem Median** gilt als nicht aktiv, denn der Wert ist nicht besser als der Median. Verglichen wird im Logarithmus mit einer Toleranz von 1e-9, damit Rundungsrauschen keine Rolle spielt.
5. **Zensierte Werte** wie in der Vorab-Spezifikation, mit dem Median als Grenze: „>" mit Grenze am oder über dem Median gilt als nicht aktiv, „>" mit Grenze unter dem Median ist nicht entscheidbar und wird ausgeschlossen. „<" mit Grenze unter dem Median gilt als aktiv, sonst ausgeschlossen (kommt nicht vor).
6. **Ausschluss.** Wird eine Verbindung nach Punkt 5 ausgeschlossen, wird der Datensatz ohne sie mit der Pipeline der Hauptanalyse neu aufgebaut (Bereinigung, Korrelationsgruppen, labelfreie Repräsentanten). Das betrifft nur DS13.
7. **Mindestgrößen.** Die Klassengrenzen der Auswahl (kleinere Klasse mindestens 25 % und 10 Verbindungen) gelten hier nicht, weil die Serien schon feststehen. Sie werden trotzdem eingehalten (siehe Tabelle).

# Folgen für die 16 Serien

| Serie | n Haupt | aktiv Haupt | Median | n Median-Regel | aktiv | kleinere Klasse | Klasse geändert | ausgeschlossen |
|---|---|---|---|---|---|---|---|---|
| DS01 | 50 | 34 | 235 nM | 50 | 25 | 50 % | 9 | 0 |
| DS02 | 34 | 14 | 50.000 nM | 34 | 14 | 41 % | 0 | 0 |
| DS03 | 45 | 12 | 5.200 nM | 45 | 22 | 49 % | 10 | 0 |
| DS04 | 40 | 23 | 869,5 nM | 40 | 20 | 50 % | 3 | 0 |
| DS05 | 42 | 30 | 2,05 µg/mL | 42 | 21 | 50 % | 9 | 0 |
| DS06 | 58 | 22 | 12,5 µg/mL | 58 | 22 | 38 % | 0 | 0 |
| DS07 | 42 | 13 | 62.500 nM | 42 | 19 | 45 % | 6 | 0 |
| DS08 | 36 | 15 | 48.000 nM | 36 | 18 | 50 % | 3 | 0 |
| DS09 | 32 | 10 | 450 nM | 32 | 16 | 50 % | 6 | 0 |
| DS10 | 45 | 19 | 213 nM | 45 | 22 | 49 % | 3 | 0 |
| DS11 | 34 | 10 | 314,5 nM | 34 | 17 | 50 % | 7 | 0 |
| DS12 | 30 | 11 | 1.000 nM | 30 | 13 | 43 % | 2 | 0 |
| DS13 | 34 | 10 | 13.345 nM | 28 | 11 | 39 % | 1 | 6 |
| DS14 | 48 | 22 | 17.420 nM | 48 | 24 | 50 % | 2 | 0 |
| DS15 | 39 | 27 | 8.900 nM | 39 | 19 | 49 % | 8 | 0 |
| DS16 | 32 | 13 | 12.450 nM | 32 | 16 | 50 % | 3 | 0 |

Lesart:

- Die Median-Regel ändert die Klasse von 0 bis 10 Verbindungen je Serie, zusammen 72 von 635 eingestuften Verbindungen.
- In DS02 und DS06 sind die Klassen identisch mit der Hauptanalyse. In DS02 und DS12 liegen viele Werte auf der höchsten getesteten Konzentration (DS02 20-mal „> 50.000 nM", DS12 16-mal „> 1.000 nM"), und der Median fällt genau auf diese Grenze. In DS06 haben 8 Verbindungen genau den Medianwert 12,5 µg/mL. Diese Verbindungen gelten als nicht aktiv, deshalb liegt der Anteil in DS02, DS06 und DS12 unter 50 %.
- DS13: 20 der 34 Werte sind zensiert. Sechs Grenzen liegen unter dem Median und sind nicht entscheidbar. Ohne sie bleiben 28 Verbindungen (1.361 Deskriptoren, 58 Gruppen).
- DS07 und DS08 sind MIC-Serien mit Werten in nM, wie in ChEMBL erfasst.

# Rechnung

- Alle zehn Analysen wie in der Hauptanalyse, mit denselben Permutationszahlen (genestete Modelle 2.000, scheinbare Modelle und Trefferzahl 1.000, Random Forest 500) und 20.000 Bootstrap-Ziehungen.
- Seeds: numpy default_rng([93, Datensatznummer, Analysecode + 300, Permutationsindex]). + 100 nach der Vorab-Spezifikation, + 200 nach der Zusatzspezifikation. Bootstrap 308.
- **DS02 und DS06** haben dieselben Klassen und Daten wie in der Hauptanalyse. Ihre Ergebnisse werden übernommen, nicht neu gerechnet, wie bei DS01, DS03, DS04 und DS15 in der Hauptanalyse. Bootstrap mit Seed 308.
- Neu gerechnet werden 14 Serien. Die vier Teilmengen der Größenvariation werden nicht gerechnet, denn die Vorab-Spezifikation sieht die Median-Regel nur für die 16 Datensätze vor.
- Skripte: datasets_m.py (Klassen, Neuaufbau DS13), run_benchmark_m.py (identisch mit run_benchmark_k.py bis auf datasets_m/, ergebnisse_m/, Analysecode + 300).

# Auswertung

- Nur beschreibend. Kein Ergebnis dieser Analyse geht in die Entscheidungen zu H1 bis H4 ein.
- Berichtet werden dieselben Kennzahlen und Tests wie für H1 bis H4 (Wilcoxon, Page, Spearman, Zählung mit Bootstrap), ohne Entscheidung und ohne Holm-Korrektur, dazu die Robustheitsprüfungen R1, R2 und AUC.
- Gegenüberstellung je Serie: Hauptanalyse gegen Median-Regel (Optimismus, Inflationsfaktoren, genestete Q² und AUC).
- Ergebnis-Arbeitsmappe im gleichen Layout wie die Hauptanalyse.

# Laufzeit

Hochrechnung aus den gemessenen Paketzeiten der bisherigen Rechnungen für dieselben Serien: etwa 15,5 Stunden bei unveränderten Laufzeiten. Eine Messung mit permutierten Klassen der Median-Regel (drei Permutationen je Route in DS01, DS07 und DS10) ergab für die genesteten Routen das 1,5- bis 3-Fache der früher gemessenen Zeit je Permutation, wohl weil bei ausgeglichenen Klassen mehr Deskriptoren den Filter passieren. Erwartet werden deshalb **etwa 18 bis 27 Stunden**, also 10 bis 15 Blöcke zu 110 Minuten. Start nur nach Freigabe durch EW.

Überwachung: kein Sicherheitsnetz, das selbst rechnet. Eine Wiederaufnahme nach einem Abbruch startet einfach den nächsten Block, denn das Skript überspringt fertige Pakete.

# Abweichungen

(leer – wird nach dem Festschreiben geführt)
