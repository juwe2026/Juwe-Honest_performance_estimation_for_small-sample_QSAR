# Revision computations (October 2026)

This folder contains the scripts that produced every number added or recomputed for the revision, and
`results/revision/` contains their outputs. `CHANGELOG_revision.md` in the archive root describes the two code
corrections (sections 1 and 7) and what changed as a result.

All scripts were run in one working directory (`/home/claude/rev`); paths are absolute and have to be adjusted when
the scripts are run elsewhere. `common.py` loads the cleaned Route A data and the 182 correlation groups and defines
the shared helpers; `pval_tie.py` defines the permutation p-value with ties counted (tolerance 1e-9) and stores every
null distribution it is given. The corrected Route A modules are those in `code/`.

Four paths in the scripts point at folders of that working directory that are not separate folders here. `/home/claude/rev/QSAR_selection_leakage` is this archive as submitted (version 1.0.0); the scripts read the Route A matrix, the Route B descriptor file and the VIP workbook from its `data/` folder; its `data/RouteB_descriptor_data.xlsx` and `data/20260710_VIP_EffectSize_Tabellen.xlsx` are `data/RouteB_descriptor_data_submitted.xlsx` and `data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx` here. `/home/claude/rev/fixcode` held the corrected modules: every module the revision scripts import from it (`nested_generic.py`, `altorder.py`, `unsup_rep.py`, `a2vip.py`) is byte-identical to the one in `code/`. `/home/claude/rev/auditcode` and `/home/claude/rev/auditcode_g` held copies of the Route A and Gaucher modules in which the fallback branch counts how often it is reached; they are used only by `audit_fallback.py` and `audit_gaucher.py` and are bundled here as `revision/auditcode/` and `revision/auditcode_g/`. `common.py` loads `base.pkl` (the five strain data sets, the 182 correlation groups, their adjacency matrix and the label-free representatives); it is not bundled because of its size (about 175 MB), and `make_base.py` rebuilds it from the code in `code/` and Additional file 6 (verified identical to the file used, object for object).

| Script | Output (`results/revision/`) | Where it is reported |
|---|---|---|
| `jobI2_fallback.py` (seeds 31, 47, 7) and `jobL_a2hi.py` (seed 42, Route A2 for *H. influenzae*) | `I_fallback_seeds.json`, `A2_HI_seed42.json` | Table 1, Routes A1, A1u, A2 |
| `jobJ2_topk.py` (seeds 59, 137, 91, 11, 23) | `cache_topk_seeds/` | Table 3, Table 4, Additional files 9 and 10, top-6 in the Discussion |
| `jobA_nofilter.py` (seed 8101) | `A_nofilter.json` | Table 1, arms without selection |
| `jobD_rf.py`, `jobD2_rf_morgan.py` (seed 8501, 500 permutations) | `D_rf.json`, `D_rf_morgan.json` | Table 1, Table S13 |
| `jobE_misc.py` (seeds 8302, 702, 8601), `jobE2_morgan_fixed.py` (seed 8602) | `E_misc.json`, `morgan2_fixed.json` | bootstrap intervals, Table S13 |
| `job_b25_cineole0.py` (seed 701, 2,000 permutations, the scheme of `code/run_routeB.py`, which it reproduces exactly for the submitted value) | `B25_cineole0_seed701.json`, `permB_25_cineole0/` | Table 2, 25 descriptors with the corrected stereocentre count of 1,8-cineole (0 instead of 2); the 25-descriptor run inside `jobE_misc.py` (seed 702) is superseded |
| `job_twostage25.py` (deterministic, algorithm of `code/lit_plsda.py`) | `twostage25_cineole0.json` | Table S16; sheets VIP_two_stage and VIP_gt1_list of Additional file 3 |
| `job_vip25_cineole0.py` (deterministic, two components) | `vip25_cineole0.json` | Figs. S3c and S4b (`figures_corrected/`), sheets VIP_Lit of Additional file 8 |
| `mk_af2.py`, `mk_af5_vip.py` | `af2_regen.json`, `af5_vip_regen.json` | Additional files 3 and 8 rebuilt from the submitted versions: RDKit 2026.03.6, stereocentre count of 1,8-cineole, correlation matrices and effect sizes with ties between isomers |
| `bootstrap_predictions_ds1.py`, `bootstrap_predictions_ds2.py` (deterministic, no permutations; random forest random_state 0) and `paired_bootstrap.py` (seed 8302, 20,000 stratified resamples per endpoint, the first 2,000 being those of `jobE_misc.py`; Holm and Benjamini-Hochberg over all 104 comparisons) | `bootstrap_predictions_ds1.json`, `bootstrap_predictions_ds2.json`, `paired_bootstrap_Q2.json`, `paired_bootstrap_Q2.csv` | Table S11 and every difference in honest Q² interpreted in Results, Discussion and Conclusions (an interval that excludes zero is required) |
| `jobC_routeB.py` (seed 8301) | `C_routeB.json` | Route B ablations (Results), Table S9 |
| `jobH_modal.py`, `jobH2_median_fixed.py`, `jobH3_median_nonnested_tie.py` (seeds 9101–9108) | `H_modal_*.json`, `H_median_fixed_*.json` | Table S3 |
| `jobF2_no3d_fixed.py` (seeds 8701, 8702) | `F_no3d_fixed.json` | response to Reviewer 2, comment 14 |
| `jobB_hits.py` (seeds 8202 and 8203; the calibrated inflation factors use no new permutations but the published null arrays of `hits_three_variants.npy`, seed 9, see `code/make_hits_three_variants.py`), `jobG_gaucher_hits.py` (seed 8801) | `B_hits.json`, `G_gaucher_hits.json` | Table S12, Table S5, Fig. 7 |
| `lone/run2_tie.py`, `lone/run3.py` (seeds 9201–9207, rule C +100) | `lone/res_*.json` | Table 5, Additional file 4 |
| `job_endpoint_corr.py` | `endpoint_corr.json` | Table S19 |
| `job_lv_noselection.py` (deterministic, no permutations) | `lv_noselection.json` | Table S8 (continued), arms without selection with one to four latent variables |
| `job_ionisation.py` (deterministic) | `ionisation_check.json` | Methods: carboxylates and tropolone tautomers (response to Reviewer 2, comment 14) |
| `job_lv_routes.py` (deterministic, no permutations; needs `QSAR_ROUTEA_XLSX` and `QSAR_ROUTEB_XLSX` pointing at the files in `data/`) | `lv_routes.json` | Table S8 (continued), Routes A1 and B with one to four latent variables: fit to all compounds (R2), apparent and nested Q2 |
| `audit_fallback.py`, `audit_gaucher.py`, `verify_morgan_fb.py` | `fallback_audit.json`, `morgan_fb_verify.json` | how often the fallback branch is reached |
| `compare_tie.py`, `cp_borderline.py` | `tie_effect.json`, `cp_borderline.json` | effect of counting ties; Table S10 |
| `check_tables_rev.py` | `numbercheck_rev.json` | every value and p-value of Tables 1 to 5 checked against its assigned source file; sheet Number_check of Additional file 2 |

`rawmic.py` extracted the replicate MIC values from the laboratory workbook, which is not part of this archive; the
extracted values are deposited in Additional file 2 (sheet MIC_replicates_31).

`results/revision_nulls/` holds every null distribution passed to `pval_tie.pval`, one file per call, with the
observed value and the p-value with and without ties counted. File names give the script, its arguments, the process
number and a running index. The runs were interrupted once by a restart of the computing environment; every arm
finished before the restart was kept from its first run, and every arm interrupted by it was repeated in full, so a
few analyses appear twice with identical values.

`fig1_v3.py` draws Fig. 1 of the revised article at print size, 174 mm wide (`cd figures && python ../revision/fig1_v3.py output`); it checks that every text lies inside its box, that no two boxes overlap and that no lettering is smaller than 8 pt.

`figures_corrected/` holds Figs. S3 and S4 as printed in the revised Supplement; the same files are now also in
`figures/output/` under the script's own file names. They are produced by `figures/figS2_figS3_vip.py` from the revised
Additional file 8, which `data/20260710_VIP_EffectSize_Tabellen.xlsx` now is, byte for byte; with the submitted workbook,
kept as `data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx`, the script reproduces the submitted figures pixel for
pixel. In the same way the Route B descriptor file in `data/` is the revised Additional file 3 and the submitted one is
kept beside it as `data/RouteB_descriptor_data_submitted.xlsx` (see the README). `20260802_FigS5_S6_Rohdaten.xlsx` is
unchanged from the submission.

`requirements-revision-lock.txt` lists the exact package versions of the environment in which the revision analyses
were run (RDKit 2026.03.6); the versions used for the submitted analyses were not recorded.
