# Changes for the second revision (November 2026)

Release 1.2.0 of the archive. Release 1.1.0 accompanied the first revision (`CHANGELOG_revision.md`), release 1.0.0 the
submitted manuscript. No analysis of releases 1.0.0 and 1.1.0 was changed, and every archived result object of those
releases is unchanged.

## 1. ChEMBL benchmark

New folder `benchmark/`: the pre-specification of 22 October 2026, the supplementary specification of 24 October 2026
and the addendum for the median rule of 25 October 2026 with their SHA-256 checksums, the code, the derived data sets
(20 each for the main analysis and the analysis under the original rule, 16 for the median rule), all permutation packages, the evaluation and robustness results, the three results workbooks as
written by the code and the logs of the runs. `benchmark/README_benchmark.md` describes the pipeline and what has been
checked. The ChEMBL-derived data are distributed under CC BY-SA 3.0 (`benchmark/DATA_LICENCE.md`). The ChEMBL database
file (2.1 GB) and the articles read for the selection are not included.

## 2. Analyses of the second revision

New folder `revision2/` with the analyses added for the second revision (Monte Carlo uncertainty of all reported
permutation p-values, Table S11, the tie-break rules against 0.0167, Table S6, the effective number of endpoints,
Table S12, the post hoc analysis of the inflation factor, Additional file 9, sections B9 and B10, the apparent Q2 of
Fig. 6a), the scripts of Fig. 1 and Fig. 6, the scripts that updated the Additional files and the SHAP checks behind
the response to the reviewers. `revision2/README_revision2.md` maps every script to its result and to the item of the
article.

## 3. Figures

- Fig. 1 is drawn by `revision2/figures/fig1_rev2c.py`, which writes `figures/output/Fig1_scheme_revised.png` and
  `.pdf`. `revision/fig1_v3.py` drew Fig. 1 of the first revision and is kept unchanged.
- Fig. 6 is new, drawn by `revision2/figures/fig6_benchmark.py` (`figures/output/Fig_ChEMBL_benchmark.png` and `.pdf`).
- The figure texts contain no semicolons. Only punctuation and capitalisation changed, the words are unchanged:
  `figures/fig2_fig4_comparison.py` (Fig. 2 and Fig. S4, reading aid, three semicolons replaced by full stops, final
  full stop added), `figures/fig6_order_comparison.py` (Fig. 4, reading aid, two semicolons replaced by full stops,
  final full stop added), `figures/fig7_figS6_hitcount.py` (Fig. 5, note a, one semicolon replaced by a full stop,
  Fig. S6 from the same script is unchanged) and `figures/figS4_figS5_scores.py` (Fig. S7 and Fig. S8, one semicolon
  each replaced by a full stop, final full stop added). The six figures were redrawn with these scripts.
- `figures/output/` holds the 14 figures of the article. Every PNG is byte-identical with the image in the manuscript or
  in the Supplement, and each is reproduced pixel for pixel by its archived script: the scripts in `figures/` with
  matplotlib 3.10.9 (`revision/requirements-revision-lock.txt`), the two scripts in `revision2/figures/` with matplotlib
  3.11.2 (`revision2/requirements-revision2-lock.txt`).

## 4. Correction of the empty-selection fallback: tables moved from Additional file 2

The four sheets that documented the correction of the fallback prediction (`CHANGELOG_revision.md`, section 1) are no
longer part of Additional file 2. They are kept, with their values before and after the correction, in
`results/revision/fallback_correction_tables.xlsx` (sheets `Fallback_correction`, `Fallback_correction_top15`,
`Fallback_correction_TableS2` and `Fallback_affected_files`). In the version of Additional file 2 of the first revision,
the seed column of the sheets `Fallback_correction` and `Fallback_correction_top15` showed 702, the seed of a first run
that was discarded. The values themselves were from the runs with the seeds 31, 47, 7 and 42 (Fallback_correction) and
59, 137, 11 and 23 (Fallback_correction_top15). The workbook gives these seeds, and all 35 rows were checked again
against the result objects. The sheet `Fallback_branch_audit` and the sheet `Seeds` remain in Additional file 2.

## 5. Notes on the first revision removed from the Additional files

The Additional files describe their content as it is and no longer how it changed. The notes on the changes of the
first revision that were removed from them are kept word for word in `revision2/history_notes_first_revision.md`.
Among them is the column `NumStereoCenters_submitted` of the Route B descriptor file, which `data/RouteB_descriptor_data.xlsx`
keeps.

## 6. Language of the code

All comments, docstrings and console messages are English. File and folder names, the keys of the result files, the
cell texts of internal working workbooks, the selection logs of the benchmark and the regular expressions that
recognise German text stay German, because the scripts read them back or write them into the archived results. The
main README explains the German terms (section "Language of the code").

The translation and the corrections of section 10 changed comments, docstrings and string literals in 147 scripts. For
every changed script the syntax
tree, with all string constants replaced by a placeholder and all expressions inside f-strings kept, was compared with
that of the script before the translation: the code is unchanged in all of them. In addition, the evaluation and
robustness scripts of the benchmark were run again and reproduce their archived results byte for byte, and all figure
scripts were run again and reproduce their figures pixel for pixel. One script depended on the wording of a comment:
`revision2/figures/fig1_rev2c.py` reads the approved wording of Fig. 1 from `revision2/figures/fig1_rev2.py` between two
comment markers, and its end marker was adapted to the translated comment.

## 7. Numbering

The article numbers its tables, figures and Additional files in the order of their first citation. The README files
use the published numbers. `CHANGELOG_revision.md` records release 1.1.0 and keeps the numbering of the first
revision. Docstrings, comments, sheet names and file names in `revision2/` use the working numbering of the revision,
which `revision2/README_revision2.md` maps to the published one. Comments and help texts of the scripts in `code/`,
`revision/` and `figures/` keep the numbering of release 1.1.0 (`CHANGELOG_revision.md`, sections 11 and 12). Only
`code/prun_a1.py` follows the published numbering and says so. File names in `figures/output/` and `results/` keep the
numbering of the submitted version, so that the scripts still run.

First revision -> published article:

| Kind | Mapping |
| --- | --- |
| main tables | 1 -> 1, 2 -> S16, 3 -> S17, 4 -> S18, 5 -> 2. Table 3 is new. |
| main figures | 1 -> 1, 2 -> 2, 3 -> 3, 4 -> S4, 5 -> S5, 6 -> 4, 7 -> 5. Fig. 6 is new. |
| supplementary tables | S1 -> S1, S2 -> S2, S3 -> S3, S4 -> S4, S5 -> S6, S6 -> S7, S7 -> S8, S8 -> S9, S9 -> S10, S10 -> S11, S11 -> S13, S12 -> S20, S13 -> S14, S14 -> S23, S15 -> S15, S16 -> S26, S17 -> S19, S18 -> S24, S19 -> S12, S20 -> S21, S21 -> S22, S22 -> S25. Table S5 is new. |
| supplementary figures | S1 -> S6, S2 -> S1, S3 -> S2, S4 -> S3, S5 -> S7, S6 -> S8 |
| Additional files | 1 -> 1, 2 -> 2, 3 -> 4, 4 -> 6, 5 -> 5, 6 -> 3, 7 -> 10, 8 -> 13, 9 -> 14, 10 -> 15, 11 -> 11, 12 -> 16, 13 -> 7. Additional files 8, 9, 12, 17 and 18 are new. |

## 8. Data files

The files in `data/` are unchanged. `data/RouteB_descriptor_data.xlsx` is Additional file 3 of the first revision and
`data/20260710_VIP_EffectSize_Tabellen.xlsx` Additional file 8 of the first revision, byte for byte. Additional files 4
and 13 of the article hold the same values. Their explanatory texts follow the published numbering, Additional file 4
no longer carries the column `NumStereoCenters_submitted` (section 5), and the sheet `FigS4c_axis_key` of Additional
file 13 is called `FigS3c_axis_key`.

## 9. Environment

`revision2/requirements-revision2-lock.txt` gives the package versions of the second revision (Python 3.13). The
version of Mordred, which computed the descriptors of the benchmark data sets, was not recorded. The descriptor matrices
are archived, and no later step needs Mordred.

## 10. Smaller corrections

- `README.md`, `figures/README.md`, `revision/README_revision.md` and `CHANGELOG_revision.md`: the 104 semicolons taken
  over from release 1.1.0 replaced by commas or full stops. Content unchanged.
- `revision2/figures/fig6_benchmark.py`: the docstring names the inputs in the published numbering.
- Usage lines in the docstrings of `revision2/analysis/` give the paths of this archive.
- `benchmark/maxz.py`: the docstring records the check of 6 November 2026.
- `results/revision/fallback_correction_tables.xlsx`, sheet Info, and the script that writes it name the script at its
  place in this archive (`revision2/tables/update_af2_rev2.py`).

## 11. Changes of 7 November 2026

After an external review of the second revision the authors approved the following changes.

- Fig. 1 (`revision2/figures/fig1_rev2.py`, `fig1_rev2c.py`): panel d reads "37 all-D peptides", because six of the 37
  peptides of the second data set have five to seven residues. The box "Same analysis plan" names the separate
  implementation of the benchmark (Methods, ChEMBL benchmark). The note above the dashed frame of panel a names the
  selection steps.
  `fig1_rev2c.py` draws the figure 186.5 mm high instead of 199.5 mm, with smaller gaps and the note above the dashed
  frame in one line, so that figure and legend fit into the 225 mm of the journal. Wording and font sizes are
  unchanged.
- Fig. 3 (`figures/fig3_fig5_figS1_permutation.py`): height 185 mm instead of 221 mm. Fig. S5 and Fig. S1 are unchanged.
- Fig. 4 (`figures/fig6_order_comparison.py`): the note reads "P. aeruginosa (three actives, descriptive only)" and is
  set in 8 pt.
- Fig. 6 (`revision2/figures/fig6_benchmark.py`): "Peptides" instead of "Octapeptides", also in the data-set names it
  reads from Additional file 2.
  The figure is 172.8 mm high instead of 193.1 mm, with panels a and b flatter, a smaller gap above panels c and d and
  the reading aid in one line ("Additional file 9, B6"), so that figure and legend fit into the 225 mm of the journal.
  Data, font sizes and line widths are unchanged.
- `figures/output/` holds the four new images. Each PNG is byte-identical with the image in the manuscript.
- Additional files 12 and 17: two cells each were edited after the translation by `translate_workbooks.py`, the data
  source in the sheet README (ChEMBL 37, database version of 1 May 2026, released on 29 May 2026) and item 7 of the sheet
  Exploratory (the median rule was computed, see Additional files 9 and 18). All other cells are unchanged.
- Additional file 2: the data-set names "Octapeptides" read "Peptides", and "(artefact)" reads "(descriptive)" in the
  tables that reproduce printed values.
- `benchmark/README_benchmark.md`: the seeds of the three analyses and the series taken over between them are named.
- `benchmark/README_benchmark.md` and `benchmark/DATA_LICENCE.md`: ChEMBL 37 is dated as in its release notes, database
  version of 1 May 2026, released on 29 May 2026.
- `revision2/analysis/mordred_check.py` with `revision2/results/Mordred_Pruefung 06.11.2026 2224.json`: recomputed
  with Mordred 2.0.7 and RDKit 2026.09.1, 55 of the 56 archived descriptor matrices of the benchmark agree within a
  relative 1e-12. DS06 of the original rule differs by up to 0.2 % in five descriptors. `mordred_ds06_check.py` and
  `mordred_ds06_models.py` give the effect on its correlation groups and models (`benchmark/README_benchmark.md`).
- `code/prun_a1.py`: mode `v1` computes the top 15 by VIP with one component in the VIP model and seed 137, the model
  of Table S18. Mode `v`, which `README.md` had named for Table S18, uses two components and seed 73 and does not
  correspond to it. `README.md` gives the commands for `v1`.
- `revision2/tables/check_af2.py`: reads the column `Descriptive_only` and the names "Peptides" and "(descriptive)" of
  Additional file 2 (13,429 checks, no deviation).
- `benchmark/README_benchmark.md`: the cleaning of the descriptors and the use of the stored correlation groups are
  described, and step 1 names source, size, checksum and command for the ChEMBL dump.
- `revision2/analysis/duplicate_nzv_audit.py` and `duplicate_groups_check.py` count the exact duplicate and
  near-zero-variance columns of the archived matrices and give the effect of the duplicates on the correlation groups.

## 12. Changes of 13 November 2026
- Paths: every script finds the archive root from its own location and runs from any working directory and wherever
  the archive is unpacked. The paths of the working environment in 93 scripts of `code/`, `benchmark/`, `revision/`
  and `revision2/` were replaced by archive-relative ones (`README.md`, Conventions). Every changed script was started
  from a copy of the archive in another place and from an unrelated working directory. Seventeen of them ran to the
  end and were compared with their archived results. Fifteen agree exactly or to within 4.4e-16, and `code/nlvcv.py`
  and `code/vip_fig.py` agree to within 3e-7 on the revised Route B file, as documented. Three faults found in this
  way were corrected. The scripts in `revision/` read the Route B descriptor file as submitted by default, as they did
  when they were run (`revision/common.py`). `revision/compare_tie.py` stops if its backup input is missing, instead
  of overwriting `tie_effect.json` with an incomplete comparison. `benchmark/kontrollen/baue_tabelle.py` opens the
  ChEMBL subset read-only, and `revision2/analysis/instabile_vorhersagen_tabelle.py` stops if it is missing, so that
  neither creates an empty database file.
- Gaucher data: `data/gaucher_source/` holds the file `gaucherdata.zip` of Smit et al. (2007) with its licence, and
  `code/gaucher_make_csv.py` rebuilds both Gaucher input files from it byte for byte. `README.md` names the source and
  the licence of these data and the spectrum left out. `revision2/analysis/gaucher_spectrum_check.py` identifies that
  spectrum by label-free criteria, and `gaucher_selection_sensitivity.py` shows that the nested results do not depend
  on which spectrum is left out.
- Tables S17 and S18: `results/raw_permutations/a1_top15_corrected/` holds the permutation draws of the top-15 models
  as the corrected code writes them. `revision2/analysis/a1_top15_draws_check.py` shows that they reproduce the
  published summaries exactly and that the older archived draws differ from them only where a fold reaches the
  fallback branch.
- `revision/make_morgan_bits.py` rebuilds the Morgan fingerprint matrix `results/revision/morgan2_bits.npy`, which was
  not archived, from the SMILES of Additional file 3. With it the Morgan analyses in `revision/` run from the archive
  and reproduce their archived results (`morgan_fb_verify.json` in full, the observed values of `D_rf_morgan.json` and
  `morgan2_fixed.json`).
- `revision2/analysis/groups_csv_roundtrip.py`: read with `float_precision='round_trip'`, the CSV files of all 56
  benchmark matrices give exactly the stored correlation groups.
- `revision2/analysis/duplicate_hits_check.py`: hit counts and inflation factors of the 16 series without the exact
  duplicate columns (descriptive, Additional file 9, section B6).
- `check_all.py` in the archive root recomputes the Gaucher input files, the hit-count null distributions, the
  observed honest Q2 and AUC of Table 1 and the three benchmark evaluations in one call, compares them with the
  archive and, given Additional files 1 and 2 and the main text, checks Additional file 2 against its sources. In the
  environment of `revision2/requirements-revision2-lock.txt` all six checks pass (Additional file 2: 15,123 checks,
  no deviation).

## 13. Further changes of 13 November 2026
- Univariate hit counts: the permutation streams of all hit-count tests of the three data sets and of the two fixed
  tie-break rules were extended from 1,000 to 100,000 permutations, with the first 1,000 permutations unchanged and
  checked (`revision2/analysis/hits_extended.py`). Tables S11 and S20 and Additional file 2 report these p-values in
  addition. Null means, standard deviations and inflation factors remain those of the first 1,000 permutations.
- Classification metrics of the control arms of Table 1 and Table S14 (`revision2/analysis/classification_controls.py`),
  added to Table S19 and Additional file 11. The out-of-bag predictions of the random forests were refitted with the
  settings of `revision/jobD_rf.py` and reproduce the archived out-of-bag Q2 and AUC exactly.
- Additional file 13: the Route B VIP columns held values in 9 rows of descriptors outside the 12-descriptor modelling
  set. These cells are now empty, and the rows are joined by descriptor name (`revision2/tables/fix_af13_routeB.py`).
- `results/revision/fallback_correction_tables.xlsx`: column Script_entry_point and a note on the numbering.
- Additional file 2: Tables S4, S11, S19 and S20 follow Additional file 1, and HitCounts_all_data_sets has seven new
  columns (`revision2/tables/update_af2_r2.py`).
- Table 3 of the main text shows the four hypotheses as columns. The sheet Main_Table_3 of Additional file 2 keeps one
  row per hypothesis, and `revision2/tables/check_af2.py` compares the two with rows and columns exchanged.
- Descriptive checks of the deposited data (`revision2/analysis/data_description_checks.py`, result in
  `revision2/results/Datenbeschreibung_Pruefung 13.11.2026 2348.json`): missing values of the Route A matrix (the three-dimensional
  PaDEL block of nine compounds, 3,195 cells filled with the column median), replicate support of the activity classes
  under the mean and the median rule, identical descriptor rows in Route B and in DS10 and DS15, right-censored values
  of the benchmark series and the molecular weight of the peptides. The main text (Methods, legend of Table 2), the
  legend of Table S3 and Additional file 9 (section B6) now state these facts. No result changes.
- Numerical form of the descriptor pool of the peptides (`revision2/analysis/input_precision_check.py`,
  `revision2/analysis/rounding_sensitivity.py`): the correlation grouping and the observed hit counts depend on the
  last digits of the descriptor values (README, "Known limitations"). The sheet Descriptors_1275 of Additional file 6
  now holds the matrix used for all results bit for bit. Before, 10,996 of its 47,175 cells differed from it in the
  last digit (relative difference up to 5.3e-16), which gave 106 instead of 105 repeated columns and other Route A2
  hit counts for three strains. No reported result changes.
- `revision2/analysis/mordred_peptides_check.py`: independent recomputation of the descriptor pool of the peptides
  from the SMILES with mordredcommunity 2.0.7 and RDKit 2026.09.1 (same 1,275 columns, largest relative difference
  1.2e-14, result in `revision2/results/Mordred_Peptide_Pruefung 13.11.2026 2348.json`). README, the lock file of the second
  revision and the Availability section of the main text state this.
- `requirements.txt`: the upper bound of pandas is raised to below 4.0, and the note refers to the two lock files.
- `revision2/analysis/review_sensitivity_checks.py` (result in
  `revision2/results/Sensitivitaet_SD-Regel_RF-Seeds_RouteB-LGO 14.11.2026 0520.json`): the rule that discards
  descriptors with a training standard deviation below 1e-6 is applied only in the arms without filter. Removing the
  two affected descriptors of the monoterpenoids and the one of the peptides changes no observed Q² or AUC of Routes A1,
  A1u and A2. The random forest of the peptides was repeated with 20 seeds, and Route B was repeated with identical
  descriptor vectors left out together. The main text states the rule and the leave-group-out result of Route B (Methods), and the legend of Table S14 gives the
  seed range. The deviations after the pre-specification of the second data set are listed in Additional file 6 (sheet
  Info).
- `revision/lone/PRESPEC.md`: the note of 1 October 2026 that fixed the activity definition of the second data set
  before any model was fitted (German, 3,266 bytes). The file itself had not been archived. It was restored from the
  unchanged German text reproduced in Additional file 5, and its SHA-256 (8b563093...51d2) equals the value recorded
  there. `check_all.py` verifies this checksum as its sixth check.
- Fig. 1 (`revision2/figures/fig1_rev2c.py`): the dashed frame no longer encloses the box of Route A1u, because the
  centrality representatives are computed once on all compounds without the labels. Only the outline of the frame
  changed. `figures/output/Fig1_scheme_revised.png` and `.pdf` are redrawn.
- Review report of 15 November 2026 (items G1, G2, G7, G10, G12 and G13, approved by EW on 16 November 2026):
  - `revision2/analysis/threshold_one_step_down.py` (result in `revision2/results/Schwelle_eine_Stufe_tiefer
    16.11.2026 1245.txt`): with the threshold one dilution step lower (mean MIC below 1024 µg/mL), carvacrol, thymol,
    thymoquinone and β-thujaplicin are active for S. aureus, S. pneumoniae and S. pyogenes, and these four with citral
    for H. influenzae. The main text (Methods) states this.
  - `revision2/analysis/rand_index_partitions.py` (result in `revision2/results/Rand_Index_Partitionen_RouteA
    16.11.2026 1245.json`): the all-singletons partition has a Rand index of 0.928 against the 182 groups of the
    monoterpenoids, because 93 % of the descriptor pairs fall between groups. The adjusted Rand index of the 31
    fold-internal partitions against the global one is 0.49 to 0.96 (median 0.90) and below 0.6 in 4 folds.
    Additional file 1 (section before Table S2) states this. Table S2 is unchanged.
  - `revision2/analysis/rounding_nested_q2.py` (result in `revision2/results/Rundung_geschachtelte_Q2_Peptide
    16.11.2026 1245.json`): the nested models of the peptides refitted on the descriptor matrix rounded to 12, 10
    and 8 significant digits, without permutations. The nested Q² of Routes A1, A1u and A2 changes by up to 0.19,
    PLS-DA on all 1,275 descriptors does not change. The main text (Methods of the second data set) and the README
    ("Known limitations") state this. Table 2 is unchanged and holds for the matrix as archived.
  - Fig. S8 (`figures/figS4_figS5_scores.py`): the line printed under the figure calls the dashed outlines 95 % data
    ellipses (normal theory) instead of confidence ellipses. They are covariance ellipses of the scores of each class
    and describe the spread of the data. `figures/output/FigS5_PLSDA_score_plots.png` and `.pdf` are redrawn with
    matplotlib 3.10.9. Only the line of text changed, and the PNG is byte-identical with the image in Additional file 1.
  - `revision2/tables/aenderungen_cg2_workbooks.py`: in Additional files 4, 14 and 15, P. aeruginosa is called
    "descriptive only" instead of an "artefact" (six text cells and eight table cells). In Additional files 12 and
    17 the 106 and 216 cells of the column Value(s) in the sheet Unstable_predictions (ChEMBL relation and value, for
    example "=7000") are stored as text instead of as formulas. No number changes. The workbooks in
    `benchmark/workbooks/` and the files of the first revision in `data/` and `results/` stay as they were written.
  - `benchmark/ergebnismappe.py`: a text beginning with "=" is written as text, so that a new run does not store the
    column Value(s) as formulas again.
- Fig. 1 (`revision2/figures/fig1_rev2c.py`, with the approved wording in `fig1_rev2.py`, approved by EW on 17 November
  2026): the box of the univariate filter reads "the |effect size| > 0.3" instead of "the |r_rb| > 0.3", whose subscript
  was printed as a literal underscore. Only this line of text changed. `figures/output/Fig1_scheme_revised.png` and
  `.pdf` are redrawn with matplotlib 3.11.2.
- Fig. 6 (`revision2/figures/fig6_benchmark.py`, approved by EW on 18 November 2026): the key names the open markers
  "comparison data sets (open)" instead of "data sets of this study (open)", as the legend of Fig. 6 does, and panel c
  has a linear y axis (0 to 15) instead of a logarithmic one. The x axis of panel c stays logarithmic. Only the key and
  panel c changed. `figures/output/Fig_ChEMBL_benchmark.png` and `.pdf` are redrawn with matplotlib 3.11.2.
- README, section "Licence" (approved by EW on 18 November 2026): new paragraph on the values taken from Table 1 of
  Lone et al. (2019), which are distributed under CC BY 4.0, not under the MIT licence. The sentence on the licence at
  the top of the README now names "the third-party data" instead of only the Gaucher data.
