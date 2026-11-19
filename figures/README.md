# Figures

Every figure of the article and its Supplementary Information is generated from the
archived result objects in `../results/` and `../benchmark/`. Nothing is drawn by hand.

```bash
cd figures
export QSAR_RAW=../results/raw_permutations
export QSAR_XLS=../data
python ../revision2/figures/fig1_rev2c.py output  # Fig. 1
python ../revision/fig1_v3.py          output   # Fig. 1 of the first revision (superseded, same file name)
python fig1_scheme.py                  output   # Fig. 1 as submitted (superseded, kept unchanged)
python fig2_fig4_comparison.py         output   # Fig. 2 and Fig. S4
python fig3_fig5_figS1_permutation.py  output   # Fig. 3, Fig. S5, Fig. S1
python fig6_order_comparison.py        output   # Fig. 4
python fig7_figS6_hitcount.py          output   # Fig. 5 and Fig. S6
python figS2_figS3_vip.py              output   # Fig. S2 and Fig. S3
python figS4_figS5_scores.py           output   # Fig. S7 and Fig. S8
```

Fig. 6 is drawn by `../revision2/figures/fig6_benchmark.py` from Additional files 2, 7 and 12 and two result files
of `../revision2/results/` (usage in `../revision2/README_revision2.md`).

The scripts in this folder reproduce their figures pixel for pixel with matplotlib 3.10.9
(`../revision/requirements-revision-lock.txt`), the two scripts in `../revision2/figures/` with matplotlib 3.11.2
(`../revision2/requirements-revision2-lock.txt`).

`output/` holds the 14 figures as printed in the article. Every PNG is byte-identical with the image in the
manuscript or in the Supplement. File names keep the numbering of the submitted version, so that the scripts
still run. The right-hand column gives the number in the article.

| File in `output/` | Figure in the article |
| --- | --- |
| `Fig1_scheme_revised.png`, `.pdf` | Fig. 1 (from `../revision2/figures/fig1_rev2c.py`, PNG at 600 dpi) |
| `Fig2_three_way_comparison` | Fig. 2 |
| `Fig3_nested_permutation_RouteA1_all5` | Fig. 3 |
| `Fig4_three_way_comparison_top15` | Fig. S4 |
| `Fig5_nested_permutation_RouteA1_top15` | Fig. S5 |
| `Fig6_order_comparison` | Fig. 4 |
| `Fig7_hitcount_inflation` | Fig. 5 |
| `FigS1_nested_permutation_HI` | Fig. S1 |
| `FigS2_VIP_per_descriptor` | Fig. S2 (VIP values of Additional file 13, read from `../data/20260710_VIP_EffectSize_Tabellen.xlsx`) |
| `FigS3_VIP_curves` | Fig. S3 (the same source) |
| `FigS4_LV1_class_separation` | Fig. S7 |
| `FigS5_PLSDA_score_plots` | Fig. S8 |
| `FigS6_hitcount_null_distributions` | Fig. S6 |
| `Fig_ChEMBL_benchmark.png`, `.pdf` | Fig. 6 (from `../revision2/figures/fig6_benchmark.py`) |

Each script writes a PNG at 300 dpi and a PDF (vector) of the same figure.

`figstyle.py` holds the shared style. The colour convention is identical in every
figure:

| | colour | |
| --- | --- | --- |
| Route A1 | `#E69F00` | orange, supervised choice of the group representative |
| Route A1u | `#009E73` | green, label-free choice of the group representative |
| Route A2 | `#0072B2` | blue, filter first, grouping second |
| Route B | `#CC79A7` | purple, a priori literature set |
| apparent estimates | pale tint of the route colour | |
| nested estimates | solid route colour | |
| null distributions | `#D9D9D9` grey | |
| observed value | solid line in the route colour | |
| 95th percentile of the null | dashed black line | |
| P. aeruginosa | hatched | artefact, three actives |
| active / non-active compounds | `#D55E00` / `#56B4E9` | Fig. S7 and Fig. S8 |

The palette is the Okabe-Ito colourblind-safe set. Adjacent-pair separation was
checked for protanopia, deuteranopia and tritanopia. Every pair used side by side
passes a ΔE of 8 or more in all three simulations. Red and green are deliberately
never used as a contrasting pair.

`dataio.py` maps each permutation file to its route:

| prefix | model |
| --- | --- |
| `perm_alt_*` | Route A1, full descriptor set |
| `perm_*` | Route A2, full descriptor set |
| `perm_uns_*` | Route A1u |
| `perm_a1e_*` | Route A1, top 15 by effect size |
| `perm15_*` | Route A2, top 15 by effect size |
| `perm_a1vip1_*` | Route A1, top 15 by VIP, one latent variable |
| `permB_13_*`, `permB_25_*` | Route B |
| `perm_a1k6_*` | Route A1, top 6 |
| `gperm*_*` | Gaucher comparison data set |


## Addendum, 2026-08-20: Route B comprises 12 descriptors

`C_Count` was removed from the Route B modelling set because it is perfectly correlated
with `EsterLacton_flag` (r = rho = 1.000). The a priori set of 25 descriptors is
unchanged. Affected: Fig. 2, Fig. 4, Fig. S2 panel b, Fig. S3 panel a and Fig. S4.

The recomputed Route B results live as JSON next to the scripts:

| File | Contents |
| --- | --- |
| `data_recalc/routeB12.json` | Q2, AUC and permutation p at nLV = 1 and nLV = 3, 2,000 permutations each |
| `data_recalc/vip_routeB.json` | Apparent VIP and the mean over the LOO folds, two latent variables |

They are loaded via `dataio.routeB12()` and `dataio.vip_routeB()`.

**Fixed along the way:** panel c of Fig. S2 and panel b of Fig. S3 previously showed only
the 13 reduced descriptors, even though the legend announced the complete a priori set of
25. The cause was a filter in `figS2_figS3_vip.py` that discarded rows without a Route B
VIP value. The panels now read the 25-descriptor model's columns and show all 25.

`../data/20260710_VIP_EffectSize_Tabellen.xlsx` is Additional file 8 of the first revision, byte for byte. Additional
file 13 of the article holds the same values, with its explanatory texts in the published numbering and the sheet
`FigS4c_axis_key` renamed `FigS3c_axis_key`. `figS2_figS3_vip.py` draws Figs. S2 and S3 as printed in the Supplement. Their 25-descriptor panels (Fig. S2c, Fig. S3b) use the
VIP values recomputed with the corrected stereocentre count of 1,8-cineole. With
`QSAR_XLS` pointing at a folder in which `../data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx` carries the name
`20260710_VIP_EffectSize_Tabellen.xlsx`, the script reproduces the submitted Figs. S2 and S3 pixel for pixel.
