# Vorab-Festlegung (vor Kenntnis der Ergebnisse) – zweiter MHK-Datensatz

Datum: 01.10.2026, vor Beginn der Rechnungen. Festgelegt von EW, umgesetzt von Claude.

## Datensatz
Lone A, Thomsen TT, Nielsen JE, Thulstrup PW, Klitgaard RN, Løbner-Olesen A, Lund R, Jenssen H,
Hansen PR. Structure-Activity Study of an All-d Antimicrobial Octapeptide D2D.
Molecules 2019;24(24):4571. doi:10.3390/molecules24244571 (CC BY).
37 Peptide (Leitpeptid D2D + 36 Analoga), vier Stämme, MHK per Bouillon-Mikrodilution in MHB-II,
Triplikate, laut Tabellenfußnote nach CLSI. Höchste getestete Konzentration 128 µg/mL;
">128" ist rechtszensiert.

## Aktivitätsdefinition
Es gibt in der Publikation KEINEN Grenzwert für aktiv/nicht-aktiv, und es existieren für diese
Peptide weder CLSI-/EUCAST-Breakpoints noch ECOFF-Werte (Breakpoints sind substanz- und
speziesspezifisch und setzen klinische Daten voraus). Der Schnittpunkt ist daher ein
datengetriebener Screening-Schnittpunkt ohne klinische Interpretation – dieselbe Situation
wie bei den Monoterpenoiden (Gutachter 2, Punkt 4).

- HAUPTANALYSE (Regel A): aktiv = MHK <= 8 µg/mL; für P. aeruginosa <= 16 µg/mL.
  Begründung: an der Verteilung je Stamm gewählt, damit beide Klassen besetzt sind;
  P. aeruginosa eine Verdünnungsstufe höher, weil die Serie dort insgesamt schwächer wirkt.
  Erwartete Klassengrößen: EC 20/37, SA 23/37, PA 14/37, AB 14/37.
- SENSITIVITÄTSANALYSE (Regel C): aktiv = MHK <= 16 µg/mL für alle vier Stämme.
  Erwartete Klassengrößen: EC 24/37, SA 26/37, PA 14/37, AB 24/37.

## Deskriptoren
Mordred (1613, nur 2D) + RDKit (217) = 1830 roh. Reinigung analog zum Manuskript:
nicht-numerisch entfernt, Spalten mit >50 % fehlend entfernt, Rest median-imputiert,
konstante Spalten entfernt -> 1275 Deskriptoren. Verhältnis n:p = 37:1275 (~1:34).
ChemDes war nicht verfügbar; das wird im Manuskript offengelegt.

## Zu rechnende Analysen (wie im Manuskript, Route B entfällt mangels Literaturdeskriptoren)
E1  Route A1  (Gruppierung zuerst, Repräsentant nach Effektgröße, dann univariater Filter), genestet
E1  Route A1u (Repräsentant nach Zentralität, labelfrei), genestet
E1  Route A2  (Filter zuerst, dann Gruppierung), genestet
E2  Ohne jede Auswahl: PLS-DA auf allen 1275 Deskriptoren
E2  Ohne Filter: PLS-DA auf den labelfreien Gruppenrepräsentanten
E4  Random Forest (300 Bäume, min_samples_leaf = 1) ohne Auswahl, OOB
E6  Nullverteilung der univariaten Trefferzahl für A2, A1, A1u; Inflationsfaktoren
    nominal (5 %) und auf die eigene Null-Trefferrate kalibriert

## Festgelegte Parameter
PLS-DA mit 1 latenter Variable; Leave-one-out; Autoskalierung nur auf den Trainingsverbindungen
jedes Folds; Spalten mit Trainings-SD < 1e-6 im Fold verworfen (labelfrei).
Univariater Filter: Mann-Whitney p <= 0.05 und |r_rb| > 0.3, in jedem Fold neu.
Korrelationsgruppierung: |Pearson| oder |Spearman| > 0.7, most-connected-first,
Gleichstand -> größte Varianz (wie im Manuskript, jetzt dokumentiert).
Permutationen: 2000 je Modell (Random Forest 500), Trefferzahl 1000.
Seeds: A1 9201, A1u 9202, A2 9203, ohne Auswahl 9204, ohne Filter 9205, RF 9206, Trefferzahl 9207;
jeweils numpy default_rng([seed, Permutationsindex]); Sensitivitätsanalyse Regel C: Seeds + 100.
