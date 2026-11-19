# Revision computations (October 2026)

This folder contains the scripts that produced every number added or recomputed for the revision, and
`results/revision/` contains their outputs. `CHANGELOG_revision.md` in the archive root describes the two code
corrections (sections 1 and 7) and what changed as a result. References to tables, figures and Additional files in this
file follow the published numbering of the article (version 1.2.0 of the archive).

All scripts were run in one working directory. In this archive they find the archive root from their own location and
run from any working directory. They read their inputs from the archive and write their results to
`results/revision/`, or to `QSAR_REV_OUT` if that is set. Intermediate files, among them `base.pkl`, go to
`results/recalc/` (or `QSAR_RECALC_DIR`). `common.py` loads the cleaned Route A data and the 182 correlation groups and defines
the shared helpers. The module `pval_tie.py` defines the permutation p-value with ties counted (tolerance 1e-9) and stores every
null distribution it is given. The corrected Route A modules are those in `code/`.

In the working directory, four folders held what is now part of this archive. The archive as submitted (version 1.0.0) supplied the Route A matrix, the Route B descriptor file and the VIP workbook. Its `data/RouteB_descriptor_data.xlsx` and `data/20260710_VIP_EffectSize_Tabellen.xlsx` are `data/RouteB_descriptor_data_submitted.xlsx` and `data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx` here, and the scripts read these two files by default (`QSAR_ROUTEB_XLSX` and, for `mk_af5_vip.py`, `QSAR_AF5_SUBMITTED` override them). A folder of corrected modules supplied `nested_generic.py`, `altorder.py`, `unsup_rep.py` and `a2vip.py`, each byte-identical to the one in `code/`, from which the scripts now import them. Two folders held copies of the Route A and Gaucher modules in which the fallback branch counts how often it is reached. They are used only by `audit_fallback.py` and `audit_gaucher.py` and are bundled here as `revision/auditcode/` and `revision/auditcode_g/`. `check_seeds.py` was run with the modules of the archive as submitted, which still had the former fallback. It now imports those of `code/` and confirms the same seeds, with a deviation of 0 in the first three null values of each of the four files it checks. `common.py` loads `base.pkl` (the five strain data sets, the 182 correlation groups, their adjacency matrix and the label-free representatives). It is not bundled because of its size (about 175 MB), and `make_base.py` rebuilds it from the code in `code/` and Additional file 3 (verified identical to the file used, object for object). A rerun stores the null distributions in `nulls/` of its output folder. Those of the archived runs are in `results/revision_nulls/`, named by script, arguments, process number and call.

Some scripts need files that are not part of this archive and stop with an error that names the missing file. `rawmic.py` reads the workbook of the replicate MIC readings (`QSAR_RAW_MIC_XLSX`, its values are Additional file 2, sheet MIC_replicates_31), and `job_endpoint_corr.py` reads its output `raw_mic_31.json`. `mk_af2.py` and `check_tables_rev.py` read the table extract `table1_pub_exact.json`, and `check_tables_rev.py` also the manuscript of the first revision. `compare_tie.py` reads the results before the treatment of ties (`strict_backup_20261003/`) and writes nothing without them. The results of all these scripts are archived in `results/revision/`.

`jobD_rf.py`, `jobD2_rf_morgan.py`, `jobE_misc.py`, `jobE2_morgan_fixed.py` and `verify_morgan_fb.py` read the Morgan fingerprint matrix `results/revision/morgan2_bits.npy` (radius 2, 2,048 bits, of which the 295 that vary across the 31 compounds are kept). The matrix used for the revision was not archived. `make_morgan_bits.py` rebuilds it from the SMILES column of Additional file 3, and the rebuilt matrix reproduces `morgan_fb_verify.json` exactly (`verify_morgan_fb.py 2000`) and the observed values of `D_rf_morgan.json` and `morgan2_fixed.json`.

| Script | Output (`results/revision/`) | Where it is reported |
|---|---|---|
| `jobI2_fallback.py` (seeds 31, 47, 7) and `jobL_a2hi.py` (seed 42, Route A2 for *H. influenzae*) | `I_fallback_seeds.json`, `A2_HI_seed42.json` | Table 1, Routes A1, A1u, A2 |
| `jobJ2_topk.py` (seeds 59, 137, 91, 11, 23) | `cache_topk_seeds/` | Table S17, Table S18, Additional files 14 and 15, top-6 in the Discussion |
| `jobA_nofilter.py` (seed 8101) | `A_nofilter.json` | Table 1, arms without selection |
| `jobD_rf.py`, `jobD2_rf_morgan.py` (seed 8501, 500 permutations) | `D_rf.json`, `D_rf_morgan.json` | Table 1, Table S14 |
| `jobE_misc.py` (seeds 8302, 702, 8601), `jobE2_morgan_fixed.py` (seed 8602) | `E_misc.json`, `morgan2_fixed.json` | bootstrap intervals, Table S14 |
| `job_b25_cineole0.py` (seed 701, 2,000 permutations, the scheme of `code/run_routeB.py`, which it reproduces exactly for the submitted value) | `B25_cineole0_seed701.json`, `permB_25_cineole0/` | Table S16, 25 descriptors with the corrected stereocentre count of 1,8-cineole (0 instead of 2). The 25-descriptor run inside `jobE_misc.py` (seed 702) is superseded |
| `job_twostage25.py` (deterministic, algorithm of `code/lit_plsda.py`) | `twostage25_cineole0.json` | Table S26, sheets VIP_two_stage and VIP_gt1_list of Additional file 4 |
| `job_vip25_cineole0.py` (deterministic, two components) | `vip25_cineole0.json` | Figs. S2c and S3b (`figures_corrected/`), sheets VIP_Lit of Additional file 13 |
| `mk_af2.py`, `mk_af5_vip.py` | `af2_regen.json`, `af5_vip_regen.json` | Additional files 4 and 13 rebuilt from the submitted versions: RDKit 2026.03.6, stereocentre count of 1,8-cineole, correlation matrices and effect sizes with ties between isomers |
| `bootstrap_predictions_ds1.py`, `bootstrap_predictions_ds2.py` (deterministic, no permutations, random forest random_state 0) and `paired_bootstrap.py` (seed 8302, 20,000 stratified resamples per endpoint, the first 2,000 being those of `jobE_misc.py`, Holm and Benjamini-Hochberg over all 104 comparisons) | `bootstrap_predictions_ds1.json`, `bootstrap_predictions_ds2.json`, `paired_bootstrap_Q2.json`, `paired_bootstrap_Q2.csv` | Table S13 and every difference in honest Q² interpreted in Results, Discussion and Conclusions (an interval that excludes zero is required) |
| `jobC_routeB.py` (seed 8301) | `C_routeB.json` | Route B ablations (Results), Table S10 |
| `jobH_modal.py`, `jobH2_median_fixed.py`, `jobH3_median_nonnested_tie.py` (seeds 9101–9108) | `H_modal_*.json`, `H_median_fixed_*.json` | Table S3 |
| `jobF2_no3d_fixed.py` (seeds 8701, 8702) | `F_no3d_fixed.json` | response to Reviewer 2, comment 14 |
| `jobB_hits.py` (seeds 8202 and 8203, the calibrated inflation factors use no new permutations but the published null arrays of `hits_three_variants.npy`, seed 9, see `code/make_hits_three_variants.py`), `jobG_gaucher_hits.py` (seed 8801) | `B_hits.json`, `G_gaucher_hits.json` | Table S20, Table S6, Fig. 5 |
| `lone/run2_tie.py`, `lone/run3.py` (seeds 9201–9207, rule C +100) | `lone/res_*.json` | Table 2, Additional file 6 |
| `job_endpoint_corr.py` | `endpoint_corr.json` | Table S12 |
| `job_lv_noselection.py` (deterministic, no permutations) | `lv_noselection.json` | Table S9 (continued), arms without selection with one to four latent variables |
| `job_ionisation.py` (deterministic) | `ionisation_check.json` | Methods: carboxylates and tropolone tautomers (response to Reviewer 2, comment 14) |
| `job_lv_routes.py` (deterministic, no permutations) | `lv_routes.json` | Table S9 (continued), Routes A1 and B with one to four latent variables: fit to all compounds (R2), apparent and nested Q2 |
| `audit_fallback.py`, `audit_gaucher.py`, `verify_morgan_fb.py` | `fallback_audit.json`, `morgan_fb_verify.json` | how often the fallback branch is reached |
| `make_morgan_bits.py` (deterministic) | `morgan2_bits.npy` | input of the Morgan fingerprint analyses, Table 1 (random forest), Table S14 |
| `compare_tie.py`, `cp_borderline.py` | `tie_effect.json`, `cp_borderline.json` | effect of counting ties, Table S11 |
| `check_tables_rev.py` | `numbercheck_rev.json` | every value and p-value of Tables 1 and 2 and Tables S16 to S18 checked against its assigned source file, sheet Number_check of Additional file 2 |

`rawmic.py` extracted the replicate MIC values from the laboratory workbook, which is not part of this archive. The
extracted values are deposited in Additional file 2 (sheet MIC_replicates_31).

`results/revision_nulls/` holds every null distribution passed to `pval_tie.pval`, one file per call, with the
observed value and the p-value with and without ties counted. File names give the script, its arguments, the process
number and a running index. The runs were interrupted once by a restart of the computing environment. Every arm
finished before the restart was kept from its first run, and every arm interrupted by it was repeated in full, so a
few analyses appear twice with identical values.

`fig1_v3.py` draws Fig. 1 of the first revision (superseded by `revision2/figures/fig1_rev2c.py`) at print size, 174 mm wide (`cd figures && python ../revision/fig1_v3.py output`). It checks that every text lies inside its box, that no two boxes overlap and that no lettering is smaller than 8 pt.

`figures_corrected/` holds Figs. S2 and S3 as printed in the Supplement. The same files are now also in
`figures/output/` under the script's own file names. They are produced by `figures/figS2_figS3_vip.py` from
`data/20260710_VIP_EffectSize_Tabellen.xlsx`, which is Additional file 8 of the first revision, byte for byte, and holds the
same values as Additional file 13 of the article. With the submitted workbook,
kept as `data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx`, the script reproduces the submitted figures pixel for
pixel. In the same way the Route B descriptor file in `data/` is Additional file 3 of the first revision (the descriptor
values of Additional file 4 of the article) and the submitted one is
kept beside it as `data/RouteB_descriptor_data_submitted.xlsx` (see the README). `20260802_FigS5_S6_Rohdaten.xlsx` is
unchanged from the submission.

`requirements-revision-lock.txt` lists the exact package versions of the environment in which the revision analyses
were run (RDKit 2026.03.6). The versions used for the submitted analyses were not recorded.
