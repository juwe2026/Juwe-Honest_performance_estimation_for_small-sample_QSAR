# -*- coding: utf-8 -*-
"""Single source of the proposed AF2 text changes (D2, D3, D4, E10). Used for the decision document and, after
approval, by update_af2_rev2.py. Every entry: (sheet, cell, old, new). 'old' is asserted before replacing."""

RN = 'results/revision/'
NU = 'results/revision_nulls/'
RP = 'results/raw_permutations/'

SEEDS = [
    # rows 2 to 5: Table 1, Routes A1, A1u, A2
    ('Seeds', 'B2', 'prun_alt.py / altorder.alt_nested_cv', 'jobI2_fallback.py / altorder.alt_nested_cv'),
    ('Seeds', 'F2', None, 'result object ' + RN + 'I_fallback_seeds.json, null distributions ' + NU + 'jobI2_fallback_x_*.npz'),
    ('Seeds', 'B3', 'prun_uns.py / unsup_rep.unsup_nested_cv', 'jobI2_fallback.py / unsup_rep.unsup_nested_cv'),
    ('Seeds', 'F3', None, 'result object ' + RN + 'I_fallback_seeds.json, null distributions ' + NU + 'jobI2_fallback_x_*.npz'),
    ('Seeds', 'B4', 'prun.py / nested_generic.nested_cv', 'jobI2_fallback.py / nested_generic.nested_cv'),
    ('Seeds', 'F4', None, 'result object ' + RN + 'I_fallback_seeds.json, null distributions ' + NU + 'jobI2_fallback_x_*.npz'),
    ('Seeds', 'B5', 'perm_runner.py / nested_h.nested_cv', 'jobL_a2hi.py / nested_generic.nested_cv'),
    ('Seeds', 'F5', None, 'result object ' + RN + 'A2_HI_seed42.json, null distributions ' + NU + 'jobL_a2hi_x_*.npz. This one strain has its own seed.'),
    # rows 6 to 10: fixed model size
    ('Seeds', 'A6', 'Route A1 top-15 by effect size (Table 3)', 'Route A1 top-15 by effect size (Table S25)'),
    ('Seeds', 'B6', 'prun_a1.py, mode e', 'jobJ2_topk.py / altorder.alt_nested_cv, 15 descriptors'),
    ('Seeds', 'F6', None, 'result objects ' + RN + 'cache_topk_seeds/*_A1_top15.json, null distributions ' + NU + 'jobJ2_topk_x_*.npz'),
    ('Seeds', 'A7', 'Route A1 top-15 by VIP, VIP from one component (Table 4)', 'Route A1 top-15 by VIP, VIP from one component (Table S26)'),
    ('Seeds', 'B7', 'prun_vip1.py', 'jobJ2_topk.py / altorder.alt_viptop, one component'),
    ('Seeds', 'F7', None, 'result objects ' + RN + 'cache_topk_seeds/*_A1_viptop1.json, null distributions ' + NU + 'jobJ2_topk_x_*.npz'),
    ('Seeds', 'B8', 'prun15.py', 'jobJ2_topk.py / nested_generic.nested_cv_a2_topk, 15 descriptors'),
    ('Seeds', 'E8', '1949-2000 as deposited, 2000 after the recomputation', 2000),
    ('Seeds', 'F8', None, 'result objects ' + RN + 'cache_topk_seeds/*_A2_top15.json, null distributions ' + NU + 'jobJ2_topk_x_*.npz'),
    ('Seeds', 'B9', 'prun_vip.py / reconstructed a2vip.py', 'jobJ2_topk.py / a2vip.nested_cv_a2_viptop, two components'),
    ('Seeds', 'E9', '1900-1982 as deposited, 2000 after the recomputation', 2000),
    ('Seeds', 'F9', None, 'result objects ' + RN + 'cache_topk_seeds/*_A2_viptop.json, null distributions ' + NU + 'jobJ2_topk_x_*.npz'),
    ('Seeds', 'B10', 'prun_a1_k.py', 'jobJ2_topk.py / altorder.alt_nested_cv, 6 descriptors'),
    ('Seeds', 'F10', 'archived file: perm_a1k6_*.npz. from the script', 'result objects ' + RN + 'cache_topk_seeds/*_A1_top6.json, null distributions ' + NU + 'jobJ2_topk_x_*.npz'),
    # rows 11, 12: arms not reported
    ('Seeds', 'F11', 'archived file: perm_a1v_*.npz. from the script; this arm is not reported', RP + 'perm_a1v_*.npz. This arm is not reported.'),
    ('Seeds', 'F12', 'archived file: perm_a1uvip_*.npz. from the script; this arm is not reported', RP + 'perm_a1uvip_*.npz. This arm is not reported.'),
    # row 13: Table S2
    ('Seeds', 'F13', 'archived file: results/recalc/fold_grouping_routeA.json. its own permutation draw, independent of Table 1, which is why individual p-values differ from that table within Monte Carlo error; the recomputation reproduces the published values exactly',
     'result object results/fold_grouping_routeA.json. Its own permutation draw, independent of Table 1, which is why individual p-values differ from that table within Monte Carlo error.'),
    # rows 14, 15: Route B
    ('Seeds', 'A14', 'Route B, 25 a priori descriptors', 'Route B, 25 a priori descriptors, stereocentre count of 2 for 1,8-cineole from the newer routine (comparison in the legend of Table S24)'),
    ('Seeds', 'F14', 'archived file: permB_25_*.npz (submitted stereocentre count of 1,8-cineole); Table 2 now uses the rerun with the corrected value below', RP + 'permB_25_*.npz'),
    ('Seeds', 'F15', 'archived file: permB_12_*.npz / permB_13_*.npz / permB_12nlv3_*.npz. from the script', RP + 'permB_12_*.npz, permB_13_*.npz, permB_12nlv3_*.npz'),
    # rows 16, 17: Gaucher
    ('Seeds', 'F16', 'archived file: gperm_*.npz. from the script', RP + 'gperm_*.npz'),
    ('Seeds', 'F17', 'archived file: gperm2_*.npz. from the script', RP + 'gperm2_*.npz'),
    # rows 18 to 21: hit counts and tie-break
    ('Seeds', 'A18', 'Hit-count null distribution, published (Table S12, Fig. 7, Fig. S1)', 'Hit-count null distribution (Table S12, Fig. 7, Fig. S1)'),
    ('Seeds', 'F18', None, RP + 'hits_three_variants.npy, written by code/make_hits_three_variants.py'),
    ('Seeds', 'C19', 'none: computed from the published null arrays (hits_three_variants.npy, seed 9 above)', 'none: computed from the null arrays of hits_three_variants.npy (seed 9, row above)'),
    ('Seeds', 'F19', 'revision; no new permutations', 'result object ' + RN + 'B_hits.json. No new permutations.'),
    ('Seeds', 'A20', 'Tie-break rule: published and first-index rule, each with its own null distribution (Table S5)', 'Tie-break rule: largest-variance and first-index rule, each with its own null distribution (Table S5)'),
    ('Seeds', 'F20', 'revision; j indexes the permutation', 'result object ' + RN + 'B_hits.json. j indexes the permutation.'),
    ('Seeds', 'C21', 'default_rng(8203), one generator drawing the 100 random rules in sequence; null distributions default_rng([8202, j])', 'default_rng(8203), one generator drawing the 100 random rules in sequence, null distributions default_rng([8202, j])'),
    ('Seeds', 'F21', 'revision; every rule is evaluated on the same 1,000 permutations', 'result object ' + RN + 'B_hits.json. Every rule is evaluated on the same 1,000 permutations.'),
    # rows 22 to 32: B gets the script, F the result object
    ('Seeds', 'B22', 'A_nofilter', 'jobA_nofilter.py'),
    ('Seeds', 'F22', 'revision', 'result object ' + RN + 'A_nofilter.json, null distributions ' + NU + 'jobA_nofilter_x_*.npz'),
    ('Seeds', 'B23', 'C_routeB', 'jobC_routeB.py'),
    ('Seeds', 'F23', 'revision', 'result object ' + RN + 'C_routeB.json, null distributions ' + NU + 'jobC_routeB_x_*.npz'),
    ('Seeds', 'A24', 'Route B, 25 descriptors, corrected stereocentre count of 1,8-cineole (Table 2)', 'Route B, 25 a priori descriptors, CIP-based stereocentre count of 1,8-cineole (Table S24)'),
    ('Seeds', 'B24', 'job_b25_cineole0.py / B25_cineole0_seed701', 'job_b25_cineole0.py'),
    ('Seeds', 'F24', 'revision; same permutation scheme as the archived permB_25_*.npz, which it reproduces exactly with the submitted value',
     'result object ' + RN + 'B25_cineole0_seed701.json, null distributions ' + RN + 'permB_25_cineole0/. Same permutation scheme as permB_25_*.npz, which it reproduces exactly with a stereocentre count of 2.'),
    ('Seeds', 'B25', 'C_routeB.boot / E_misc.nested_obs', 'jobC_routeB.py / jobE_misc.py'),
    ('Seeds', 'F25', 'revision; stratified resampling of the (class, prediction) pairs', 'result objects ' + RN + 'C_routeB.json and E_misc.json. Stratified resampling of the (class, prediction) pairs.'),
    ('Seeds', 'B26', 'D_rf', 'jobD_rf.py / jobD2_rf_morgan.py'),
    ('Seeds', 'F26', 'revision; 500 permutations for run time, p resolution 0.002', 'result objects ' + RN + 'D_rf.json and D_rf_morgan.json, null distributions ' + NU + 'jobD_rf_x_*.npz and jobD2_rf_morgan_x_*.npz. 500 permutations for run time, p resolution 0.002. The forest itself uses random_state 0.'),
    ('Seeds', 'B27', 'E_misc.morgan2', 'jobE_misc.py / jobE2_morgan_fixed.py'),
    ('Seeds', 'F27', 'revision', 'result objects ' + RN + 'E_misc.json and morgan2_fixed.json, null distributions ' + NU + 'jobE_misc_x_*.npz and jobE2_morgan_fixed_x_*.npz'),
    ('Seeds', 'A28', 'Routes without three-dimensional descriptors', 'Routes without three-dimensional descriptors (Table S23)'),
    ('Seeds', 'B28', 'F_no3d', 'jobF2_no3d_fixed.py'),
    ('Seeds', 'F28', 'revision; A1 and A1u respectively', 'result object ' + RN + 'F_no3d_fixed.json, null distributions ' + NU + 'jobF2_no3d_fixed_x_*.npz. Seeds for A1 and A1u respectively.'),
    ('Seeds', 'B29', 'G_gaucher_hits', 'jobG_gaucher_hits.py'),
    ('Seeds', 'F29', 'revision', 'result object ' + RN + 'G_gaucher_hits.json'),
    ('Seeds', 'A30', 'Median activity rule, sensitivity analysis', 'Median activity rule, sensitivity analysis (Table S3)'),
    ('Seeds', 'B30', 'H_modal_*', 'jobH_modal.py / jobH2_median_fixed.py / jobH3_median_nonnested_tie.py'),
    ('Seeds', 'F30', 'revision', 'result objects ' + RN + 'H_modal_*.json and H_median_fixed_*.json, null distributions ' + NU + 'jobH2_median_fixed_*.npz and jobH3_median_nonnested_tie_*.npz'),
    ('Seeds', 'B31', 'lone/res_A.json; Routes A1, A1u, A2 from lone/res_A_fix.json (corrected fallback)', 'lone/run2_tie.py / lone/run3.py'),
    ('Seeds', 'F31', 'revision; 500 for the random forest, 1,000 for the hit count', 'result objects ' + RN + 'lone/res_A.json, Routes A1, A1u and A2 from ' + RN + 'lone/res_A_fix.json, null distributions ' + NU + 'run2_tie_A_*.npz and run3_A_*.npz. 500 permutations for the random forest, 1,000 for the hit count.'),
    ('Seeds', 'B32', 'lone/res_C.json; Routes A1, A1u, A2 from lone/res_C_fix.json (corrected fallback)', 'lone/run2_tie.py / lone/run3.py'),
    ('Seeds', 'F32', 'revision', 'result objects ' + RN + 'lone/res_C.json, Routes A1, A1u and A2 from ' + RN + 'lone/res_C_fix.json, null distributions ' + NU + 'run2_tie_C_*.npz and run3_C_*.npz. 500 permutations for the random forest, 1,000 for the hit count.'),   # D8 (approved 08:43)
]
# new row 33 (paired bootstrap of Table S11, not listed so far)
SEEDS_NEW_ROW = ['Paired bootstrap of differences in Q² (Table S11)', 'paired_bootstrap.py', 'default_rng(8302)', '8302',
                 '20,000 resamples per endpoint',
                 'result objects ' + RN + 'paired_bootstrap_Q2.json and paired_bootstrap_Q2.csv. Stratified resampling, the first 2,000 resamples are those of the bootstrap intervals (jobE_misc.py).']

INFO = [
    ('Info', 'B7', 'one latent variable unless stated; leave-one-out; autoscaling on the training compounds of each fold; columns with a training-fold SD below 1e-6 are dropped within the fold (label-free, affects 2 of 2,208 columns)',
     'one latent variable unless stated, leave-one-out, autoscaling on the training compounds of each fold. Columns with a training-fold SD below 1e-6 are dropped within the fold (label-free, affects 2 of 2,208 columns)'),
    ('Info', 'B8', 'Mann-Whitney U with p <= 0.05 and rank-biserial |r| > 0.3; the effect-size threshold never binds, so the filter is effectively a p-value filter (see response to Reviewer 2, comment 8)',
     'Mann-Whitney U with p <= 0.05 and rank-biserial |r| > 0.3. The effect-size threshold never binds, so the filter is effectively a p-value filter'),
    ('Info', 'B10', 'P. aeruginosa has three actives of 31 and gives degenerate null distributions; its values are reported but not interpreted',
     'P. aeruginosa has three actives of 31 and gives degenerate null distributions. Its values are reported but not interpreted'),
    ('Info', 'B11', 'one row per individual MIC reading (triplicates in independent experiments; 6 to 12 values per compound, strain and phase, 9 for most); 2048 codes no inhibition at the highest concentration tested (1024 ug/mL; 256 ug/mL for thymoquinone) and is right-censored; the class threshold MIC < 2048 is applied to the per-compound value used in the models',
     'one row per individual MIC reading (triplicates in independent experiments, 6 to 12 values per compound, strain and phase, 9 for most). 2048 codes no inhibition at the highest concentration tested (1024 ug/mL, 256 ug/mL for thymoquinone) and is right-censored. The class threshold MIC < 2048 is applied to the per-compound value used in the models'),
    ('Info', 'B14', 'Lone et al., Molecules 2019;24:4571 (CC BY); activity definition pre-specified before any model was fitted',
     'Lone et al., Molecules 2019, 24, 4571 (CC BY). Activity definition pre-specified before any model was fitted'),
]

# cells repeated in a column: (sheet, column letter, old, new) applied to every cell with exactly 'old'
REPEATED = [
    ('SecondDataSet_models', 'A', 'main (MIC <= 8; <= 16 for P. aeruginosa)', 'main (MIC <= 8, <= 16 for P. aeruginosa)'),
    ('Table1_new_all_arms', 'K', 'only three actives; degenerate null distribution, p not interpretable', 'only three actives, degenerate null distribution, p not interpretable'),
    ('Endpoint_correlation', 'B', 'effective number of endpoints: Li and Ji 3.00; Nyholt/Cheverud 3.41', 'effective number of endpoints: Li and Ji 3.00, Nyholt/Cheverud 3.41'),
    ('Endpoint_correlation', 'B', 'effective number of endpoints: Li and Ji 3.00; Nyholt/Cheverud 2.73', 'effective number of endpoints: Li and Ji 3.00, Nyholt/Cheverud 2.73'),
    ('Table1_new_all_arms', 'B', 'Route B, 25 descriptors (Table 2, stereocentres of 1,8-cineole = 0)', 'Route B, 25 descriptors (Table S24, stereocentres of 1,8-cineole = 0)'),
    ('TieBreak_sensitivity', 'A', 'variance (as published)', 'largest variance (used in this work)'),
]

# D4: column headers of Null_distribution_audit (meaning verified on all 76 rows: blow-up = |Q2| > 10)
NULLAUDIT = [('Null_distribution_audit', 'A1', 'Datei', 'File'),
             ('Null_distribution_audit', 'B1', 'n_perm', 'Permutations'),
             ('Null_distribution_audit', 'C1', 'blowups', 'Permutations_with_abs_Q2_above_10'),
             ('Null_distribution_audit', 'D1', 'anteil_blowup', 'Share_with_abs_Q2_above_10_pct'),
             ('Null_distribution_audit', 'E1', 'Q2_median', 'Q2_median'),
             ('Null_distribution_audit', 'F1', 'Q2_mittel_ohne_blowup', 'Q2_mean_without_abs_Q2_above_10'),
             ('Null_distribution_audit', 'G1', 'AUC_mittel', 'AUC_mean')]

SHEET_RENAME = [('Table1_new_all_arms', 'Table1_all_arms')]

# E10: label of the Table S5 row in Table S10 and in the Monte Carlo sheets
LABEL_S5 = ('Route A1u hit count, tie-break largest variance, as published', 'Route A1u hit count, tie-break largest variance')   # D9 (approved 08:43)
TABLE_S5_A2 = ('Largest variance, as published (182 groups): hits (p)', 'Largest variance, used in this work (182 groups): hits (p)')
