# Figures

Every figure of the article and its Supplementary Information is generated from the
archived result objects in `../results/`. Nothing is drawn by hand.

```bash
cd figures
export QSAR_RAW=../results/raw_permutations
export QSAR_XLS=../data
python fig1_scheme.py                  output   # Fig. 1
python fig2_fig4_comparison.py         output   # Fig. 2 and Fig. 4
python fig3_fig5_figS1_permutation.py  output   # Fig. 3, Fig. 5, Fig. S1
python fig6_order_comparison.py        output   # Fig. 6
python fig7_figS6_hitcount.py          output   # Fig. 7 and Fig. S6
python figS2_figS3_vip.py              output   # Fig. S2 and Fig. S3
python figS4_figS5_scores.py           output   # Fig. S4 and Fig. S5
```

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
| active / non-active compounds | `#D55E00` / `#56B4E9` | Fig. S4 and Fig. S5 |

The palette is the Okabe-Ito colourblind-safe set. Adjacent-pair separation was
checked for protanopia, deuteranopia and tritanopia; every pair used side by side
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
unchanged. Affected: Fig. 2, Fig. 4, Fig. 6, Fig. S2 panel b and Fig. S3 panel a.

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
