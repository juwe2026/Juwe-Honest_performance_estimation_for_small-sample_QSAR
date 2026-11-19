# Result map

For every table and figure of the article: how it is obtained, with which seed, and where the archived result is. Paths
are relative to the archive root. Scripts of the Route A analyses need Additional file 3 placed as described in
`README.md` ("One input has to be placed by hand"). The scripts in `revision/` write to `results/revision/`, where the
archived results are, so set `QSAR_REV_OUT` to an empty folder to compare a rerun with them. Seeds are the arguments of
`numpy.random.default_rng`, with `j` the index of the permutation. "Archived object" means that the archive holds the
result but no driver script that writes it. The routines it was computed with are named where they exist. Additional
file 2 reproduces every table cell by cell, and its sheet Seeds lists every seed with its result object.

## Main text

| Item | Script and call | Seed | Archived result |
| --- | --- | --- | --- |
| Table 1, Routes A1, A1u, A2, nested | `python code/run_nested_observed.py` (observed values), `python revision/jobI2_fallback.py` and `python revision/jobL_a2hi.py` (permutations) | A1 `[31, j]`, A1u `[47, j]`, A2 `[7, j]`, A2 for *H. influenzae* `[42, j]` | `results/revision/I_fallback_seeds.json`, `A2_HI_seed42.json` |
| Table 1, arms without selection | `python revision/jobA_nofilter.py` | `[8101, j]` | `results/revision/A_nofilter.json` |
| Table 1, random forest | `python revision/jobD_rf.py` (reads `results/revision/morgan2_bits.npy`, rebuilt by `python revision/make_morgan_bits.py`) | `[8501, j]`, 500 permutations, forest `random_state` 0 | `results/revision/D_rf.json` |
| Table 1, Route B | `python code/run_routeB12.py` | `[702, j]` | `results/raw_permutations/permB_12_*.npz`, `permB_12nlv3_*.npz`, summary `figures/data_recalc/routeB12.json` |
| Table 1, apparent estimates (block c) | archived object, routine `code/ordinary_perm.py` | | `results/raw_permutations/ordinary_perm_results.pkl` |
| Table 2 | `python revision/lone/run2_tie.py A` and `python revision/lone/run3.py A` (`C` for the sensitivity definition), descriptor pool `revision/lone/desc_pool.csv` | `[9201, j]` to `[9207, j]` (A1, A1u, A2, no selection, no filter, random forest, hit count), definition C + 100 | `results/revision/lone/res_A.json`, `res_A_fix.json` (`res_C.json`, `res_C_fix.json`) |
| Table 3 | `benchmark/`, steps 5 to 9 of `benchmark/README_benchmark.md` (`run_benchmark_k.py`, `praed_k.py`, `maxz.py`, `auswertung_k.py`, `robustheit_k.py`) | `[93, series, analysis code + 200, j]`, series DS01, DS03, DS04 and DS15 taken over with + 0 | `benchmark/ergebnisse_k/auswertung.json`, `robustheit.json` |
| Fig. 1 | `python revision2/figures/fig1_rev2c.py figures/output` | | `figures/output/Fig1_scheme_revised.png`, `.pdf` |
| Fig. 2 | `cd figures && python fig2_fig4_comparison.py output` (environment variables as in `figures/README.md`) | | `figures/output/Fig2_three_way_comparison` |
| Fig. 3 | `cd figures && python fig3_fig5_figS1_permutation.py output` | as Table 1 | `figures/output/Fig3_nested_permutation_RouteA1_all5` |
| Fig. 4 | `cd figures && python fig6_order_comparison.py output` | | `figures/output/Fig6_order_comparison` |
| Fig. 5 | `cd figures && python fig7_figS6_hitcount.py output` | as Table S20 | `figures/output/Fig7_hitcount_inflation` |
| Fig. 6 | `python revision2/figures/fig6_benchmark.py <AF12.xlsx> <AF2.xlsx> <AF7.xlsx> <scheinbar.json> <endlich_n.json> <outstem>` (`revision2/README_revision2.md`) | | `figures/output/Fig_ChEMBL_benchmark.png`, `.pdf` |

## Additional file 1 (Supplementary Information)

| Item | Script and call | Seed | Archived result |
| --- | --- | --- | --- |
| Table S1 | `python code/fold_reduction.py 2000` | `[702, j]` | `results/fold_reduction.json` |
| Table S2 | `python code/fold_grouping_routeA.py 2000` | `[1234, j]` | `results/fold_grouping_routeA.json` |
| Table S2, adjusted Rand index (text before the table) | `python revision2/analysis/rand_index_partitions.py . <work_dir> <out.json>` | | `revision2/results/Rand_Index_Partitionen_RouteA 16.11.2026 1245.json` |
| Table S3 | `python revision/jobH_modal.py <strain>`, `python revision/jobH2_median_fixed.py <strain>`, `python revision/jobH3_median_nonnested_tie.py <strain>` | `[9101, j]` to `[9108, j]` | `results/revision/H_modal_*.json`, `H_median_fixed_*.json` |
| Table S4 | no computation (comparison with the reference methods) | | |
| Table S5 | `python revision/jobF2_no3d_fixed.py` | A1 `[8701, j]`, A1u `[8702, j]` | `results/revision/F_no3d_fixed.json` |
| Table S6 | `python revision/jobB_hits.py`, `python revision2/analysis/jobB_tiebreak_0167.py <archive_root> <out.json>` | `[8202, j]`, random rules 8203 | `results/revision/B_hits.json`, `revision2/results/TableS5_TieBreak_0.0167 06.11.2026 0334.json` |
| Table S7 | `python code/desc_compute.py`, check `python code/add_structure_identifiers.py data/RouteB_descriptor_data.xlsx` | | `data/RouteB_descriptor_data.xlsx` (Additional file 4) |
| Table S8 | sheet Reduction of `data/RouteB_descriptor_data.xlsx`, read by `code/core.py` (`RED12`). `code/reduction.py` recomputes the grouping from the published correlation matrices (see "Known limitations") | | Additional file 4, sheet Reduction |
| Table S9, Route A1 | archived objects | | `results/raw_permutations/mono_nlv.pkl`, `mono_dcv.pkl` |
| Table S9, Route B | `python code/extra_tables.py`, p-value at three components `python code/run_routeB12.py` | `[702, j]` | `results/extra_tables.json`, `results/raw_permutations/routeB_nlv3.pkl`, `permB_12nlv3_*.npz` |
| Table S9 (continued) | `python revision/job_lv_routes.py`, `python revision/job_lv_noselection.py` (deterministic) | | `results/revision/lv_routes.json`, `lv_noselection.json` |
| Table S10 | `python revision/jobE_misc.py` (Routes A1, A1u, A2), `python revision/jobC_routeB.py` (Route B) | bootstrap 8302, permutations `[8301, j]` | `results/revision/E_misc.json`, `C_routeB.json` |
| Table S11 | `python revision2/analysis/cp_borderline_all_rev2.py <archive_root> <AF2.xlsx> <out.json>`, `python revision2/analysis/hits_extended.py mono 100000` and `… tie 100000` (column p (10⁵)) | | `revision2/results/TableS10_alle_pWerte 06.11.2026 0626.json`, `hits_extended_monoterpenoids 13.11.2026 1559.json`, `hits_extended_tiebreak 13.11.2026 1647.json` |
| Table S12 | `python revision/job_endpoint_corr.py` (needs the replicate MIC readings, see `revision/README_revision.md`), independently `python revision2/analysis/meff_endpunkte.py <folder of the Additional files> <out.json>` | | `results/revision/endpoint_corr.json`, `revision2/results/Effektive_Endpunkte_Nachrechnung 05.11.2026 1555.json` |
| Table S13 | `python revision/bootstrap_predictions_ds1.py`, `python revision/bootstrap_predictions_ds2.py`, `python revision/paired_bootstrap.py` | 8302, 20,000 resamples per endpoint | `results/revision/paired_bootstrap_Q2.json`, `.csv` |
| Table S14 | `python revision/jobD_rf.py`, `jobD2_rf_morgan.py`, `jobE_misc.py`, `jobE2_morgan_fixed.py` (Morgan matrix as for Table 1) | `[8501, j]` (500), `[8601, j]`, `[8602, j]` | `results/revision/D_rf.json`, `D_rf_morgan.json`, `E_misc.json`, `morgan2_fixed.json` |
| Table S15 | archived object | | `results/raw_permutations/r2_vs_q2.json` |
| Table S16 | 25 descriptors `python revision/job_b25_cineole0.py`, 12 descriptors `python code/run_routeB12.py` | `[701, j]`, `[702, j]` | `results/revision/B25_cineole0_seed701.json`, `results/raw_permutations/permB_12_*.npz` |
| Table S17 | apparent: archived object. Nested: `python code/prun_a1.py e <SHORT> Cleaned_Steps1-4 2000 <SHORT>` | `[59, j]` | `results/raw_permutations/eff15_apparent.pkl`, `results/revision/cache_topk_seeds/<SHORT>_A1_top15.json`, draws `results/raw_permutations/a1_top15_corrected/perm_a1e_<SHORT>.npz` |
| Table S18 | apparent: archived object. Nested: `python code/prun_a1.py v1 <SHORT> Cleaned_Steps1-4 2000 <SHORT>` | `[137, j]` | `results/raw_permutations/vip1lv_apparent.pkl`, `vip1lv_final.npy`, `results/revision/cache_topk_seeds/<SHORT>_A1_viptop1.json`, draws `results/raw_permutations/a1_top15_corrected/perm_a1v1_<SHORT>.npz` |
| Table S19 | compiled from the confusion matrices (TP, FN, FP, TN) that the result objects of the models carry, and for the control arms `python revision2/analysis/classification_controls.py` | | Additional file 11, `revision2/results/classification_controls 13.11.2026 1636.json` |
| Table S20 | `python code/make_hits_three_variants.py` (null hit counts), `python revision/jobB_hits.py` (calibrated factors), `python revision/jobG_gaucher_hits.py` (Gaucher data), `python revision2/analysis/hits_extended.py mono 100000` and `… gau 100000` (columns p (10⁵)) | 9, new generator per strain, Gaucher `[8801, j]` | `results/raw_permutations/hits_three_variants.npy`, `results/revision/B_hits.json`, `G_gaucher_hits.json`, `revision2/results/hits_extended_monoterpenoids 13.11.2026 1559.json`, `hits_extended_gaucher 13.11.2026 1559.json` |
| Table S21 | observed values `code/gaucher_pipe.py`, permutations `python code/gperm.py <route> 2000` and `python code/gperm2.py <route> 2 2000`, doubly nested row `python code/gaucher_doublecv.py` | `[909, j]`, `[424, j]` | `results/raw_permutations/gaucher_res.pkl`, `gperm_*.npz`, `gperm2_*.npz`, `gaucher_doublecv.pkl` |
| Table S22 | `python code/gaucher_downsample.py 2000`, `python code/gaucher_bootstrap.py 2000` | `[20260829, k, d]` | `results/raw_permutations/gaucher_downsample_2000.json`, `gaucher_bootstrap_2000.json` |
| Table S23 | archived object | | `results/raw_permutations/scaling_leakage.json` |
| Table S24 | archived objects, `code/ga_run.py` (see "Known limitations") | | `results/raw_permutations/ga_HI.pkl`, `ga_SA.pkl` |
| Table S25 | archived objects, routines in `code/lpo_cv.py` | | `results/raw_permutations/lpo_all_routes.json`, `lpo_compare.pkl` |
| Table S26 | `python revision/job_twostage25.py` (deterministic) | | `results/revision/twostage25_cineole0.json` |
| Fig. S1 | `cd figures && python fig3_fig5_figS1_permutation.py output` | as Table 1 | `figures/output/FigS1_nested_permutation_HI` |
| Figs. S2, S3 | `cd figures && python figS2_figS3_vip.py output` | | `figures/output/FigS2_VIP_per_descriptor`, `FigS3_VIP_curves` |
| Fig. S4 | `cd figures && python fig2_fig4_comparison.py output` | | `figures/output/Fig4_three_way_comparison_top15` |
| Fig. S5 | `cd figures && python fig3_fig5_figS1_permutation.py output` | as Table S17 | `figures/output/Fig5_nested_permutation_RouteA1_top15` |
| Fig. S6 | `cd figures && python fig7_figS6_hitcount.py output` | as Table S20 | `figures/output/FigS6_hitcount_null_distributions` |
| Figs. S7, S8 | `cd figures && python figS4_figS5_scores.py output` | | `figures/output/FigS4_LV1_class_separation`, `FigS5_PLSDA_score_plots` |
