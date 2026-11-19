# Second revision (November 2026)

This folder holds the scripts and results of every computation added for the second revision, apart from the ChEMBL
benchmark itself, which is in `../benchmark/`. All numbers below follow the published numbering of the article.

## Numbering in the scripts

Docstrings, comments, sheet names and file names in this folder were written during the revision and use its
working numbering: the numbering of the first revision, with Fig. 8, Tables S23 to S26 and Additional files 14 to 18
for the new items. `tables/nummern_freigegeben.py` maps every working number to the published one:

| Kind | Working number -> published number |
| --- | --- |
| main tables | 1 -> 1, 5 -> 2, 6 -> 3. The former main Tables 2, 3 and 4 -> Tables S16, S17 and S18. |
| main figures | 1 -> 1, 2 -> 2, 3 -> 3, 6 -> 4, 7 -> 5, 8 -> 6. The former main Figs. 4 and 5 -> Figs. S4 and S5. |
| supplementary tables | S1 -> S1, S2 -> S2, S3 -> S3, S4 -> S4, S5 -> S6, S6 -> S7, S7 -> S8, S8 -> S9, S9 -> S10, S10 -> S11, S11 -> S13, S12 -> S20, S13 -> S14, S14 -> S23, S15 -> S15, S16 -> S26, S17 -> S19, S18 -> S24, S19 -> S12, S20 -> S21, S21 -> S22, S22 -> S25, S23 -> S5 |
| supplementary figures | S1 -> S6, S2 -> S1, S3 -> S2, S4 -> S3, S5 -> S7, S6 -> S8 |
| Additional files | 1 -> 1, 2 -> 2, 3 -> 4, 4 -> 6, 5 -> 5, 6 -> 3, 7 -> 10, 8 -> 13, 9 -> 14, 10 -> 15, 11 -> 11, 12 -> 16, 13 -> 7, 14 -> 9, 15 -> 8, 16 -> 12, 17 -> 17, 18 -> 18 |

In usage lines, `<RouteA_matrix.xlsx>` is Additional file 3 (placed as described in `../README.md`) and
`<archive_root>` is the root folder of this archive.

## Environment

`requirements-revision2-lock.txt` (Python 3.13). Fig. 1 and Fig. 6 are reproduced pixel for pixel with matplotlib
3.11.2 from this file.

## Scripts and results

Results are in `results/` unless stated otherwise. The time stamp in a file name is the time of the run (Berlin time).

### Analyses

| Script | Result | Used for |
| --- | --- | --- |
| `analysis/cp_borderline_all_rev2.py` | `TableS10_alle_pWerte 06.11.2026 0626.json` | Table S11 (Monte Carlo uncertainty of all reported permutation p-values, Clopper-Pearson intervals, thresholds 0.05 and 0.0167), sheet MonteCarlo_uncertainty of Additional file 2 |
| `analysis/cp_borderline_rev2.py` | none kept | first version of the script above, restricted to the main-text tables, superseded |
| `analysis/jobB_tiebreak_0167.py` | `TableS5_TieBreak_0.0167 06.11.2026 0334.json` | Table S6 (the 100 random tie-break rules of the grouping against the threshold 0.0167), Additional file 2 |
| `analysis/meff_endpunkte.py` | `Effektive_Endpunkte_Nachrechnung 05.11.2026 1555.json` | Table S12 (effective number of independent endpoints), independent recomputation |
| `analysis/scheinbare_q2_rev1code.py` | `Scheinbare_Q2_Octapeptide_Monoterpenoide 05.11.2026 2059.json` | Fig. 6a (apparent minus nested Q2 of the monoterpenoids and peptides) |
| `analysis/naeherung_endlich_n.py` | `Inflation_Naeherung_endliches_n_posthoc 05.11.2026 2232.json` | Additional file 9, section B9, Table B9, Fig. 6d |
| `analysis/inflation_naeherung.py` | `Inflation_Naeherung_posthoc 05.11.2026 2216.json` | Additional file 9, section B9 (leading term of the approximation) |
| `analysis/studien_naeherung.py` | `Inflation_Naeherung_Studiendaten_posthoc 05.11.2026 2225.json` | Additional file 9, section B9, the three data sets of the study |
| `analysis/h3_sekundaer_pruefung.py` | `H3_sekundaer_Pruefung 05.11.2026 2033.json` | Additional file 9, section B9 (why G/p does not predict the factor), input of the two approximation scripts |
| `analysis/gaucher_indikatorkorrelation.py` | console output | Additional file 9, Table B8 |
| `analysis/gegenpruefung_benchmarkcode.py` | `Gegenpruefung_Benchmarkcode_Studiendaten 05.11.2026 2111.json` | Additional file 9, section B10, Table B10 (benchmark code against the code of this archive) |
| `analysis/vergleich.py` | `Vergleich_Rev1code_Benchmarkcode 05.11.2026 2111.json` | the same comparison per endpoint and route |
| `analysis/benchmark_kennzahlen_text.py` | `Benchmark_Kennzahlen_Ergebnistext 05.11.2026 2012.json` | key figures of the Results section on the benchmark that are not pre-specified (post hoc, descriptive) |
| `analysis/instabile_vorhersagen_tabelle.py` | internal working workbook (German) | the series DS03, DS11 and DS15 with the unstable predictions, compound by compound, as in Additional file 12, sheet Unstable_predictions |
| `analysis/gruppe_test.py` | console output | check that the vectorised correlation grouping gives the same groups as `../benchmark/datasets.py` |
| `analysis/mordred_check.py` | `Mordred_Pruefung 06.11.2026 2224.json` | descriptor matrices of the benchmark recomputed with Mordred 2.0.7 (the version used at the time was not recorded) |
| `analysis/mordred_ds06_check.py` | `Mordred_DS06_Pruefung 07.11.2026 1915.json` | DS06 of the original rule, the one matrix that differs: descriptors, ranks, correlation graph and groups |
| `analysis/mordred_ds06_models.py` | `Mordred_DS06_Modelle 07.11.2026 1915.json` | DS06 of the original rule: all observed values on the recomputed matrix |
| `analysis/duplicate_nzv_audit.py` | `Duplikate_NZV_Audit 13.11.2026 1106.json` | exact duplicate and near-zero-variance columns in the archived matrices (Methods, second data set and benchmark) |
| `analysis/duplicate_groups_check.py` | `Duplikate_Gruppen 13.11.2026 1106.json` | correlation groups and G/p of DS01 to DS16 with and without the duplicate columns |
| `analysis/a1_top15_draws_check.py` | `A1_top15_Ziehungen_v1 13.11.2026 1206.json`, `A1_top15_Ziehungen_e 13.11.2026 1211.json` | Tables S17 and S18: the draws of `../code/prun_a1.py` against the published summaries and the older archived draws (README.md, Known limitations) |
| `analysis/groups_csv_roundtrip.py` | `Gruppen_CSV_Rundlauf 13.11.2026 1041.json` | correlation groups of all 56 benchmark matrices regrouped from the CSV files, with the default and with the round-trip float parser of pandas |
| `analysis/duplicate_hits_check.py` | `Duplikate_Treffer 13.11.2026 1102.json` | hit counts, calibrated inflation factors and the tests of H2 and H3 of DS01 to DS16 without the exact duplicate columns (descriptive) |
| `analysis/gaucher_spectrum_check.py` | `Gaucher_Spektren 13.11.2026 1211.json` | the 40 spectra of the published Gaucher file compared by label-free criteria, to identify the one left out (../README.md, The Gaucher comparison data) |
| `analysis/gaucher_selection_sensitivity.py` | `Gaucher_Auswahl 13.11.2026 1212.json` | nested results of the Gaucher data with all 40 spectra and with each patient spectrum left out in turn (Table S21) |
| `analysis/hits_extended.py` | `hits_extended_monoterpenoids 13.11.2026 1559.json`, `hits_extended_peptides 13.11.2026 1558.json`, `hits_extended_gaucher 13.11.2026 1559.json`, `hits_extended_tiebreak 13.11.2026 1647.json` | permutation p-values of every univariate hit count with the permutation streams extended to 100,000, the first 1,000 permutations checked against the archived ones (Tables S6, S11 and S20, Additional file 2) |
| `analysis/mordred_peptides_check.py` | `Mordred_Peptide_Pruefung 13.11.2026 2348.json`, `Mordred_Peptide_je_Deskriptor 13.11.2026 2348.csv` | descriptor pool of the 37 peptides recomputed from the SMILES with Mordred and RDKit and compared with `revision/lone/desc_pool.csv` (independent run, needs Mordred) |
| `analysis/input_precision_check.py` | `Eingabegenauigkeit_Peptide 13.11.2026 2348.json` | correlation groups and observed hit counts of the peptides for the matrix as the scripts read it, read exactly, from Additional file 6 and rounded to 15 to 6 significant digits |
| `analysis/rounding_sensitivity.py` | `Rundungsempfindlichkeit 13.11.2026 2348.txt` | the same sensitivity for the monoterpenoids, and the number of groups of DS01 to DS16 after rounding to 10 significant digits |
| `analysis/rounding_nested_q2.py` | `Rundung_geschachtelte_Q2_Peptide 16.11.2026 1245.json` | nested Q² and AUC of the peptides (Routes A1, A1u, A2, PLS-DA on all 1,275 descriptors and on the 50 representatives) at full precision and rounded to 12, 10 and 8 significant digits, without permutations (main text, Methods of the second data set) |
| `analysis/threshold_one_step_down.py` | `Schwelle_eine_Stufe_tiefer 16.11.2026 1245.txt` | active monoterpenoids per strain when the threshold is moved one dilution step down (mean MIC below 1024 µg/mL), read from Additional file 2 (main text, Methods) |
| `analysis/rand_index_partitions.py` | `Rand_Index_Partitionen_RouteA 16.11.2026 1245.json` | Rand index and adjusted Rand index of the 31 fold-internal correlation partitions of the monoterpenoids against the global one, and of the all-singletons partition (Additional file 1, section before Table S2) |
| `analysis/data_description_checks.py` | `Datenbeschreibung_Pruefung 13.11.2026 2348.json` | missing values of the Route A matrix, replicate support of the activity classes, identical descriptor rows (Route B, benchmark), right-censored values, molecular weight of the peptides |
| `analysis/review_sensitivity_checks.py` | `Sensitivitaet_SD-Regel_RF-Seeds_RouteB-LGO 14.11.2026 0520.json` | observed Q² and AUC of Routes A1, A1u and A2 with and without the descriptors whose standard deviation is below 1e-6, out-of-bag Q² and AUC of the random forest of the peptides for 20 seeds, Route B with identical descriptor vectors left out together |
| `analysis/classification_controls.py` | `classification_controls 13.11.2026 1636.json` | confusion matrices of the control arms of Table 1 and Table S14, random forests from their out-of-bag predictions (Table S19, Additional file 11) |

`results/Hauptanalyse_auswertung 25.10.2026 1142.json`, `Hauptanalyse_maxz_labelfrei`, `Hauptanalyse_robustheit` and
`Hauptanalyse_vorhersagen` are byte-identical copies of the corresponding files in `../benchmark/ergebnisse_k/`, kept
under the names under which they were delivered.

### Figures

| Script | Output | Figure |
| --- | --- | --- |
| `figures/fig1_rev2c.py <output folder>` | `Fig1_scheme_revised.png` (600 dpi), `.pdf` | Fig. 1. It reads the wording of every box from `figures/fig1_rev2.py` and stops if a text differs. |
| `figures/fig6_benchmark.py <AF12.xlsx> <AF2.xlsx> <AF7.xlsx> <scheinbar.json> <endlich_n.json> <outstem>` | `<outstem>.png` (300 dpi), `.pdf` | Fig. 6 |
| `figures/fig1_rev2.py` | draft of Fig. 1 with panels c and d side by side, superseded | source of the approved wording |
| `figures/abb_pruefung.py` | `abb_pruefung.json` | check of the main figures against the figure requirements of the journal (size, line widths, font sizes) |

Both figures are in `../figures/output/` (`Fig1_scheme_revised` and `Fig_ChEMBL_benchmark`). They were checked pixel
for pixel against the figures of the article. `fig6_benchmark.py` gives the same image with the Additional files of
the article and with those of the first revision.

### Additional files

| Script | Purpose |
| --- | --- |
| `tables/aenderungen_af2.py`, `tables/update_af2_rev2.py`, `tables/af2_p2.py`, `tables/renumber_af2.py` | Additional file 2: content updates (Table S11 with all p-values, Table S6 against 0.0167, the effective endpoints), the new sheet Main_Table_3 and the published numbering |
| `tables/check_af2.py` | checks every sheet of Additional file 2 against its source: the result objects of this archive, Additional file 1 and the main text |
| `tables/aenderungen_af3_13.py`, `tables/apply_af3_13.py` | Additional files 3, 4, 6, 7, 11 and 13 to 16 from their versions of the first revision: titles, cross-references and sheet names in the published numbering. Every number cell is checked to be unchanged. |
| `tables/uebersetzung_af16_18.py`, `tables/translate_workbooks.py` | Additional files 12, 17 and 18 from the results workbooks of the benchmark (see `../benchmark/README_benchmark.md`) |
| `tables/nummern.py`, `tables/nummern_freigegeben.py` | the numbering in citation order and the substitution of references |
| `tables/fix_af13_routeB.py` | Additional file 13: Route B VIP columns joined by descriptor name, exactly 12 rows per strain |
| `tables/update_af11_controls.py` | Additional file 11: control arms of Table 1 and Table S14 added |
| `tables/update_af2_r2.py` | Additional file 2: Tables S4, S11, S19 and S20 rewritten from Additional file 1, seven new columns of HitCounts_all_data_sets |
| `tables/fallback_script_column.py` | `results/revision/fallback_correction_tables.xlsx`: column Script_entry_point |
| `tables/aenderungen_cg2_workbooks.py` | Additional files 4, 14 and 15: P. aeruginosa called "descriptive only" instead of "artefact". Additional files 12 and 17: the cells of the column Value(s) in the sheet Unstable_predictions stored as text instead of as formulas. Every other cell is checked to be unchanged |

The scripts in `tables/` take the Additional files as arguments. Additional file 2 of the first revision is not part of
this archive.

### SHAP (response to Reviewer 1)

`shap/figure_r1_shap_beeswarm.py` and `shap/figure_r1_shap_beeswarm_zwei_skalen.py` draw Figure R1 of the response to
the reviewers. `shap/shap_ms/` and `shap/shap_pruef/` hold the SHAP checks behind the response (interventional and
conditional SHAP values for PLS-DA with one latent variable). These results are not part of the article.

## Paths and inputs outside the archive

Every script finds the archive root from its own location and runs from any working directory. Intermediate files go
to `results/recalc/` of the archive root (or `QSAR_RECALC_DIR`). Two groups of scripts need inputs that are not in
this archive. `analysis/instabile_vorhersagen_tabelle.py` needs `../benchmark/chembl37_subset.db` from step 1 of the
benchmark pipeline and stops if it is missing. The SHAP scripts that use the data of this study
(`shap/figure_r1_shap_beeswarm.py`, `shap/figure_r1_shap_beeswarm_zwei_skalen.py` and those in `shap/shap_ms/`) read
the Route A matrix and the second data set from the folder `QSAR_REV1_AF_DIR`, under their names in the first
revision, `Additional_file_6_RouteA_descriptor_matrix 06.10.2026.xlsx` and `Additional_file_4_second_MIC_data_set
06.10.2026.xlsx`. Additional files 3 and 6 of the article can be placed there under these names. With them
`shap/shap_ms/shap_ms.py` reproduces `shap_ms_interventionell.json` exactly. All SHAP scripts need the package `shap`
(0.54.1 in `requirements-revision2-lock.txt`). `shap/shap_pruef/bedingt.py` samples without a fixed seed, so a rerun
agrees with `bedingt_DS01.json` only within the sampling variation. `figures/fig1_rev2c.py` reads only
`figures/fig1_rev2.py`, relative to its own location.
