# ChEMBL benchmark

This folder holds the code, the derived data sets and the complete results of the ChEMBL benchmark (manuscript,
section "Generalisation across 20 ChEMBL data sets", Fig. 6, Additional files 8, 9, 12, 17 and 18). It is the working
directory in which the benchmark was computed, kept as it was run, so that every script runs in place. All commands
below are run from this folder.

The data come from ChEMBL release 37 (database version of 1 May 2026, released on 29 May 2026). Licence and
attribution are in `DATA_LICENCE.md`.

## Specifications

The selection rules, the analysis plan and the four hypotheses were frozen before any model was fitted. The frozen
documents and their SHA-256 checksums are in `specifications/` and are published as Additional file 8.

| File | Frozen | SHA-256 |
| --- | --- | --- |
| `Vorab-Spezifikation_ChEMBL-Benchmark_FESTGESCHRIEBEN 22.10.2026 1810.md` (pre-specification) | 22 October 2026, 18:10 | `a13bb2c0b513289b160de744bd8d6bc47662d46553ab5a2fe3db6d336a0ae419` |
| `Zusatzspezifikation_Hauptanalyse_FESTGESCHRIEBEN 24.10.2026 1720.md` (supplementary specification, main analysis) | 24 October 2026, 17:20 | `36668775daafd821c1871a15c6a697ea5e3c25afebbb3d6b6a4db818c902df87` |
| `Ergaenzung_Median-Regel_FESTGESCHRIEBEN 25.10.2026 2139.md` (addendum, median rule) | 25 October 2026, 21:39 | `eaf21b7dff370396663d4fe93943a2fc4bd02044e9a45f00ffef94ef1b48eeab` |

The checksum files next to them (`PRESPEC_SHA256.txt`, `ZUSATZSPEZ_SHA256 24.10.2026 1720.txt`,
`ERGAENZUNG_MEDIAN_SHA256 25.10.2026 2139.txt`) allow the check `sha256sum -c`.

## Three analyses

| Analysis | Rule | Data sets | Results | Published as |
| --- | --- | --- | --- | --- |
| main analysis | supplementary rule: reference, control and comparison compounds of other classes and values taken from other publications are not part of the series | `datasets_k/` | `ergebnisse_k/` | Additional file 12 |
| analysis under the original rule | pre-specification of 22 October 2026 | `datasets/` | `ergebnisse/` | Additional file 17 |
| sensitivity analysis | median rule (class boundary at the median of the series), 16 series | `datasets_m/` | `ergebnisse_m/` | Additional file 18 |

Each data-set folder holds, per data set `DSxx`, the descriptor matrix `DSxx_X.csv`, the compounds with ChEMBL
identifier, SMILES, class and activity value `DSxx_y.csv`, and `DSxx_info.json` (cell of the design, document, assay,
target, unit, class sizes, number of descriptors before and after cleaning, correlation groups). DS01 to DS16 are the
natural series, DS21 to DS24 the four subsets of the size variation. Each results folder holds the observed values
(`beobachtet.json`), the permutation packages (`perm/`, one file per data set, model and block of permutations), the
observed predictions (`vorhersagen.json`), the label-free check of numerical stability (`maxz_labelfrei.json`), the
hypothesis tests (`auswertung.json`) and the robustness checks (`robustheit.json`).

The descriptors are computed and cleaned by `deskriptoren()` and `reinigen()` in `datasets.py`, as fixed in the
pre-specification: the two-dimensional Mordred descriptors and the RDKit descriptor list (1,830 columns),
non-numerical values set to missing, columns with more than 50 % missing values removed, the remaining gaps imputed
with the median, and columns with a standard deviation of at most 1e-12 and repeated column names removed. Columns
with identical values are not removed, and the near-zero-variance filter of the monoterpenoid Route A is not applied.
In DS01 to DS16, 96 to 140 columns per series duplicate another column, and each of them lies in the correlation group
of that column. Every later step reads the correlation groups stored in `DSxx_info.json` (key `gruppen`) and does not
regroup. Regrouping the CSV files with `gruppieren()` gives the same number of groups in 13 of the 16 series and 94,
84 and 92 instead of 95, 82 and 93 groups in DS01, DS02 and DS11. The differences come from the default float parser of pandas, which does not read every stored value back to the last
bit. Read with `float_precision='round_trip'`, the CSV files of all 56 matrices in `datasets/`, `datasets_k/` and
`datasets_m/` give exactly the stored groups (`../revision2/analysis/groups_csv_roundtrip.py`).

## Pipeline

| Step | Script | Input | Output |
| --- | --- | --- | --- |
| 1 | `stream_parse.py chembl37_subset.db` | ChEMBL 37 MySQL dump on standard input, see below | `chembl37_subset.db` (not included, 2.1 GB) |
| 2 | `candidates.py` | `chembl37_subset.db` | `cand_all.json`, `cand_ranked.json` |
| 3 | `select.py` (original rule), `auswahl_neu.py` (supplementary rule) | candidates, reading verdicts `verdikte.json`, decisions per compound `kontrollen/entscheidungen.json` | cell assignment, `kontrollen/auswahl_neu.json` with the selection log |
| 4 | `datasets.py` and `teilmengen.py` (original rule), `datasets_k.py` (main analysis), `datasets_m.py` (median rule) | `chembl37_subset.db`, step 3 | `datasets/`, `datasets_k/`, `datasets_m/` |
| 5 | `run_benchmark.py`, `run_benchmark_k.py`, `run_benchmark_m.py` | data sets | `beobachtet.json`, `perm/` |
| 6 | `praed.py`, `praed_k.py`, `praed_m.py` | data sets | `vorhersagen.json` |
| 7 | `maxz.py <data-set folder> <output file>` | data sets | `maxz_labelfrei.json` |
| 8 | `auswertung.py`, `auswertung_k.py`, `auswertung_m.py` | steps 5 to 7 | `auswertung.json` |
| 9 | `robustheit.py`, `robustheit_k.py`, `robustheit_m.py` | steps 5 to 8 | `robustheit.json` |
| 10 | `ergebnismappe.py haupt|vorab <file>`, `ergebnismappe_m.py <file>` | steps 3 to 9 | results workbooks in `workbooks/` |
| 11 | `../revision2/tables/translate_workbooks.py` | `workbooks/` | Additional files 12, 17 and 18 |

Step 1 reads the dump `chembl_37/chembl_37_mysql/chembl_37_mysql.dmp` of `chembl_37_mysql.tar.gz` (2,101,694,224
bytes, SHA-256 `e6eb871ee46121404d0caa79ebc00d451f5bd9b4b14688e223659cca216ec467` as listed in `checksums.txt`) from
https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/, for example with

```bash
tar -xzOf chembl_37_mysql.tar.gz chembl_37/chembl_37_mysql/chembl_37_mysql.dmp | python3 stream_parse.py chembl37_subset.db
```

The three nested routes are in `routen.py`, the vectorised univariate filter in `core_fast.py`. `core_fast.py` checks
itself against `scipy.stats.mannwhitneyu` (`pruefe_core_fast()`). `permtest.py` measured the run time with permuted
labels only, before the main run (`laufzeit.json`, `laufzeit2.json`).

Step 5 computes in packages and writes a file after every package, so that an interrupted run resumes where it
stopped: `python3 run_benchmark_k.py 110` computes for 110 minutes and can be called again. The logs of the runs are
`block*.log` (original rule), `hblock*.log` (main analysis) and `mblock*.log` (median rule). Nested models were
permuted 2,000 times, apparent models and the hit count 1,000 times and the random forest 500 times. The seeds follow
from the analysis codes of the specifications: analysis code + 0 for the original rule (`ergebnisse/`), + 200 for the
main analysis (`ergebnisse_k/`) and + 300 for the median rule (`ergebnisse_m/`). As specified, series that are
identical between two analyses were taken over and not recomputed: DS01, DS03, DS04 and DS15 in `ergebnisse_k/` come
from `ergebnisse/` (codes + 0), and DS02 and DS06 in `ergebnisse_m/` come from `ergebnisse_k/` (codes + 200). A fresh
run of `run_benchmark_k.py` or `run_benchmark_m.py` on an empty results folder uses + 200 or + 300 for every series
and therefore reproduces these six series only after their permutation packages are copied from the earlier analysis.

Step 10 writes some cell texts in German. Step 11 translates them and checks that no cell other than text changes.
`workbooks/` holds the three workbooks of step 10 as they were written: `ChEMBL_Benchmark_Ergebnisse_Hauptanalyse
25.10.2026 2106.xlsx`, `ChEMBL_Benchmark_Ergebnisse_Vorab-Analyse 25.10.2026 2106.xlsx` and
`ChEMBL_Benchmark_Ergebnisse_Median-Regel 30.10.2026 1414.xlsx`. Additional file 12 is the translated main-analysis
workbook with two cross-references renumbered to the published numbering (sheet Deviations, cell D11, and sheet
Excluded_compounds, cell A30: Additional file 14 of the working numbering is Additional file 9).

## What has been checked

On 6 November 2026, in the environment of `../revision2/requirements-revision2-lock.txt`:

- Steps 8 and 9 reproduce `auswertung.json` and `robustheit.json` of all three analyses byte for byte.
- Step 7 reproduces `maxz_labelfrei.json` of the main analysis and of the median rule byte for byte, and that of the
  original rule to within a relative difference of 4e-13 (no series changes its flag).
- Step 11 reproduces Additional files 17 and 18 cell for cell from `workbooks/`, and Additional file 12 apart from the
  two renumbered cross-references.

During the main analysis, packages computed independently in different sessions were compared byte for byte wherever
they overlapped (22,188 comparisons), without a single difference.

## What is not included, and what that means

- `chembl37_subset.db` (2.1 GB). Steps 1 to 4 and the sheet Unstable_predictions of step 10 need it. Every later step
  reads the archived data sets, so steps 5 to 9 and 11 run from this archive alone.
- The articles that were read to decide on reference and control compounds and on the authors' thresholds, for reasons
  of copyright. The decisions are recorded in `verdikte.json`, `kontrollen/entscheidungen.json` and in Additional
  file 9, Table B1, with the source passage for every excluded compound. `dateinamen.txt` lists the file names of the
  articles, `artikel/` holds the scripts that fetched and searched the open-access texts.
- The Mordred version used in step 4 was not recorded. The descriptor matrices it produced are archived. Recomputed from
  the archived structures with mordredcommunity 2.0.7 and RDKit 2026.09.1 (`../revision2/analysis/mordred_check.py`),
  all 56 matrices have the same columns, and 55 agree within a relative difference of 1e-12. DS06 of the original rule
  (`datasets/`) differs by up to 0.2 % in RD_AvgIpc, RD_Ipc, VR1_A, VR2_A and VR3_A, which changes 3 of its 66
  correlation groups. On this recomputed matrix the hit counts, the descriptors selected on all compounds and the AUC
  of Routes A1, A1u and A2 are unchanged, their nested Q2 changes by at most 0.0007, that of the no-filter control by
  0.007 (AUC 0.005), the random forest is unaffected, and PLS-DA on all descriptors remains numerically unstable
  (`../revision2/analysis/mordred_ds06_models.py`).

## Working scripts outside the pipeline

`abgleich.py`, `beschaffung.py`, `leseliste.py`, `scanpdf.py`, `build_xlsx.py`, `build_xlsx2.py`, `endbericht.py`,
`postprocess.py`, `profil.py` and `kontrollen/baue_tabelle.py` organised and searched the reading of the articles and
wrote internal working workbooks. Five of them (`beschaffung.py`, `build_xlsx.py`, `build_xlsx2.py`, `endbericht.py`, `kontrollen/baue_tabelle.py`)
write their workbooks to `../results/recalc/` (or `QSAR_RECALC_DIR`). The first four run from this archive.
`kontrollen/baue_tabelle.py` also needs `chembl37_subset.db` from step 1, which it opens read-only. None of them is
needed for a result.

## Language

File names, folder names, dictionary keys and the cell texts of the internal working workbooks are German, as they
were when the analyses ran. Translating them would change the archived result files, which the scripts reproduce byte
for byte. Comments, docstrings and console messages are English. The German terms are explained in the glossary of
the main README (`../README.md`, section "Language of the code").
