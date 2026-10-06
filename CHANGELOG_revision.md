# Changes to the deposited code for the revision (October 2026)

Release 1.1.0 of the archive (release 1.0.0 accompanied the submitted manuscript). The title of the article now ends "and its label-free remedy".

## 1. Corrected fallback prediction when no descriptor passes the filter

In the nested cross-validation routines there is a branch for folds in which the
univariate filter retains no descriptor, or in which the training fold contains
only one class. The prediction for the held-out compound was set to the mean of
the training labels:

    yhat[i] = ytr.mean()

For leave-one-out this equals (S - y_i) / (n - 1), where S is the sum of all
labels. It is therefore a strictly decreasing function of the held-out label, so
active compounds systematically receive the smaller prediction and such folds
contribute a reversed ranking. It is a label leak in the fallback branch. The
fallback is now the neutral constant:

    yhat[i] = 0.5

The line occurred in nine files, all of which have been corrected:
`nested_generic.py`, `altorder.py`, `unsup_rep.py`, `fold_grouping_routeA.py`,
`nested_h.py`, `prun_a1u_vip.py`, `gaucher_pipe.py`, `lit_plsda.py`, `lpo_cv.py`,
and in the reconstructed `a2vip.py` (see section 5).

### How often the branch is actually reached

We counted the folds that take it. On the observed labels it is reached in none
of the 31 folds, in any of the 45 arm-and-endpoint combinations examined, so no
published point estimate is affected; every one of them reproduces exactly to four decimal
places. Under permuted labels it is reached in 0.0 to 3.8 per cent of folds,
most often under Route A1u and for the three-active P. aeruginosa endpoint. The
counts per arm are in Additional file 2, sheet `Fallback_branch_audit`.

Two of the nine files are provably unaffected. In `lit_plsda.py` the branch is
dead code: the Route B models are fitted without a selector, so all columns are
used, and the two-stage VIP selector returns all columns when no descriptor
exceeds VIP = 1. In `lpo_cv.py` the branch assigned the same value to both
members of the pair, so the pair counted as a tie and contributed 0.5 to the
concordance, exactly as the neutral constant does; Table S22 is unchanged. In
`gaucher_pipe.py` the branch is never reached even under permutation (0 of 3,900
folds in 100 permutations per route), because 19 of the 39 samples are active;
Table S20 and Additional file 13 are unchanged.

### Recomputation

Every affected permutation test was recomputed with the seed and the permutation
count of the published run, so that the difference between published and
corrected value is the effect of the correction and not of a different
permutation stream. The seeds were first established by regenerating the first
null values of each archived file with the deposited code and comparing them to
the archived array; every assignment agrees exactly. They are listed in
Additional file 2, sheet `Seeds`.

The result: 52 of the 70 recomputed p-values are unchanged to four decimal
places, the largest change is 0.007 (Route A1u, P. aeruginosa, reported descriptively only), and no
verdict changes anywhere. Tables 3 and 4 reproduce exactly, p-value for p-value.
In Table S2, sixteen of the twenty p-values reproduce exactly and four change in
the third decimal place. Sheets `Fallback_correction`,
`Fallback_correction_top15` and `Fallback_correction_TableS2` of
Additional file 2 list every published and corrected value side by side, and sheet
`Fallback_affected_files` gives the account file by file.

The error was real and had to be fixed. Its effect on this particular data set is
negligible, because the affected permutations lie far below the observed value in
any case; it would not have been negligible had the observed and the null
distributions lain closer together.

Arms that contain no supervised filter (the unselected controls, Route B and the
random forest) have no such branch and are unaffected.

## 2. Added script for the hit-count null distribution

The script generating `hits_three_variants.npy` (Table S12, Fig. 7, Fig. S1) was
missing from the submitted archive. It has been reconstructed as
`code/make_hits_three_variants.py`. The permutation scheme of the archived object is a new
generator `numpy.random.default_rng(9)` for each strain and its 1,000 successive calls of
`permutation(y)`; with it the script reproduces all five observed triples and all fifteen
null arrays element for element, and the file it writes is byte-identical to the archived
one, so Table S12, Fig. 7 and Fig. S1 are unchanged. A first draft of the script used
`default_rng([8201, j])`, which reproduces the observed counts but not the null arrays; it
was never used for a reported number. The script builds its inputs (the 182 groups and the
label-free representatives) from the archive code itself and stores the Route A1 count under
the key `A1s`, as the archived object does. The inflation factors calibrated to each route's
own null hit rate (Table S12, continued) are computed from these archived null arrays by
`revision/jobB_hits.py`; that script draws new permutations only for the tie-break analysis
of Table S5 (seeds 8202 and 8203).

## 3. Tie-break rule of the correlation grouping documented

The grouping resolves ties in the most-connected-first ordering by the larger
variance. This was implemented but not stated in the manuscript; it is now
documented in the Methods and in the code comments. Resolving ties by first index
instead yields 184 groups, which accounts exactly for the count reported by
Reviewer 2.

## 4. Permutation counts

The files `perm15_*.npz` contain 1,949 to 2,000 permutations and `perm_vip_*.npz`
1,900 to 1,982, rather than the nominal 2,000, because permutations that produced
a degenerate training fold were dropped. All four top-15 arms have been rerun to
exactly 2,000 permutations; the actual count is stated with each result in
Additional file 2.

## 5. Missing function reconstructed

`prun_vip.py`, the script behind the Route A2 VIP arm of Additional file 10, calls
`nested_generic.nested_cv_viptop`. That function is not present in the deposited
module, so the arm could not be reproduced from the archive as deposited. It has
been reconstructed as `a2vip.nested_cv_a2_viptop`, line for line as
`nested_cv_a2_topk` but selecting by VIP instead of by effect size in the final
selection step. It reproduces all ten archived values (Q2 and AUC for five
strains) exactly.

The reconstruction also revealed a documentation error: the VIP ranking of this
arm is computed from a two-component PLS-DA, not from the one-component model
stated in the Info sheet of Additional file 10. The call omits the `nlv_vip`
argument and the default is two. With one component the values differ by up to
0.11 in Q2. The Route A1 arm of Table 4 does use one component, as documented
(`prun_vip1.py` passes `nlv_vip=1` explicitly). The Info sheet has been
corrected and the number of components is now stated for every VIP calculation.

`prun_vip.py` itself now calls the reconstructed function, with the same arguments
as the original call (VIP from two components). It also names the activity column
explicitly, as `prun_a1.py` does, because the pooled sheet `Cleaned_Steps1-4` carries
all five activity columns, and it writes its output archive-relative to
`results/recalc/`. Run as `python code/prun_vip.py HI Cleaned_Steps1-4 2000 HI`, it
reproduces the observed Q2 and AUC of all five strains exactly and the first 20 null
values of each strain exactly, except two values of *P. aeruginosa* (permutations 11
and 16). In each of these one training fold reaches the fallback branch, which now
predicts 0.5 (section 1). With the former fallback these two values agree exactly too.

## 6. Seeds documented

The seed of every permutation run is now listed, with its script, its archived
file and its permutation count, in Additional file 2, sheet `Seeds`. One
detail worth noting: the Route A2 run for H. influenzae was produced by a
separate script (`perm_runner.py`, seed 42) while the other four strains come
from `prun.py` (seed 7). The grouping-sensitivity analysis of Table S2 uses its
own draw (seed 1234), which is why its p-values differ from those of Table 1
within Monte Carlo error for the same quantity.

## 7. Ties in the permutation p-values

`p = (1 + #{permuted >= observed}) / (1 + n_permutations)` counts permuted values equal to the observed one as at least
as extreme. The deposited scripts compared with a strict floating-point `>=`. The AUC of a leave-one-out prediction
for n1 actives and n0 non-actives is a fraction k/(n1 n0); the same fraction computed from different predictions is
stored with differences of the order of 1e-16 (for example 102/130 as 0.7846153846153845, ...846 and ...847), so part
of the ties were not counted and some AUC p-values came out slightly too small. In the archived null distributions this
concerns 19 of 76 files; Q2 values are continuous and are not affected in practice.

All comparisons now use a tolerance of 1e-9 (`null >= observed - 1e-9`), in every script under `code/` that computes a
permutation p-value. Every reported AUC p-value of the revision has been recomputed with this rule, using the published
seeds; the revision runs additionally store their null distributions (`results/revision_nulls/`). Point estimates and
Q2 p-values are unchanged (apart from two Q2 p-values of the degenerate Morgan nested-filter arm for S. pyogenes, where
many folds predict the same constant). The largest change in the main-text tables is 0.0075 (Route A1u, S. aureus, AUC,
0.0505 to 0.0580); the largest overall is 0.0115 (Table S3, H. influenzae under the median rule, Route A1, AUC). Two
nominal verdicts change, both for values within Monte Carlo error of 0.05: Additional file 10, Route A2 top-15 by effect
size, S. aureus, AUC 0.0497 to 0.0525, and Table S3, S. pyogenes under the median rule, Route A1u, AUC 0.0465 to 0.0520.
The significance markers of Fig. 2 and Fig. 4 are drawn from `figures/data_recalc/revision_pvalues.json`; one marker of
Fig. 4 (S. aureus, Route A2 nested AUC) disappears.

`results/results.json` is the output of `code/run_all.py` as submitted and was computed with the strict comparison.
Re-running the script now gives the tie-corrected values. This was checked from a clean unpacking of this archive
(4 October 2026): 676 of 680 entries agree exactly; the four that differ are two AUC permutation p-values and their
counts. Route B with 12 descriptors and one latent variable, S. pneumoniae: 0.0065 to 0.0075 (12 to 14 permuted values
at least as large); Table 1 reports 0.0075. The superseded 13-descriptor set with three latent variables, S. pyogenes:
0.1114 to 0.1124; this set is not reported in the article.

## 8. Figures
- Fig. 7a now shows the inflation of the null standard deviation relative to the binomial value both at the route's own
  null hit rate (bars) and at the nominal 5 % rate (ticks) (`figures/fig7_figS6_hitcount.py`).
- Fig. 2 and Fig. 4 take their significance markers from the recomputed p-values (see section 7).
- Fig. 1 was redrawn for the revision by `revision/fig1_v3.py` (output `figures/output/Fig1_scheme_revised.png` and `.pdf`). `figures/fig1_scheme.py`, which drew the submitted Fig. 1, is kept unchanged; its output `Fig1_two_route_scheme` was removed from `figures/output/` so that this folder holds only the figures of the revised article.
- Lettering only, data unchanged: the rotated strain names of Fig. 6 and Fig. 7 are smaller (7.5 and 8 pt instead of 10.5 pt), as are their y-axis tick labels (8 pt); the legend of Fig. 7 is split into one legend per panel, the two observed-count markers stacked above panel b; in Fig. S5 the y-axis label and tick labels, in Fig. S6 the axis labels, tick labels and the class legend are smaller (`figures/fig6_order_comparison.py`, `figures/fig7_figS6_hitcount.py`, `figures/figS4_figS5_scores.py`).
- `figures/output/` previously still held the submitted Fig. 4 and Fig. 7. Both were regenerated with the archived scripts; Fig. 1 to Fig. 7 in `figures/output/` are now identical pixel for pixel with the figures of the revised manuscript.

- Fig. 1: the correlation grouping of the survivors under Route A2 is now drawn with a green fill and an orange outline and reads “input selected with labels”. The step itself uses no labels, but its input is the output of the label-dependent filter, which is why it lies inside the validation loop while the grouping of Routes A1 and A1u lies outside it. The legend gained a corresponding swatch. Fig. 1 is drawn at print size for the Springer artwork guidelines: 174 mm wide, all lettering at least 8 pt (Arial-metric font); for that the box texts are broken over more lines, E5 to E7 are corner marks of their input boxes like E1 to E4, and the colour key has two rows. Content, colours and arrows are unchanged, except that in panel d the Morgan2 arm with the nested univariate filter is now drawn as a separate orange branch with its own PLS-DA box, since it is label-dependent and nested and therefore a sensitivity analysis rather than a control arm (`revision/fig1_v3.py`).
- Fig. 6 drawn on the same canvas as Fig. 7 (6.3 x 3.35 in, same margins and panel spacing) and saved with the same horizontal crop (`XSPAN_FIG6_FIG7` in `figures/figstyle.py`), so that at the same printed width the y-axes of both figures lie on the same vertical lines and their x-axes have the same length; legend, axis labels, tick labels and notes in the font sizes of Fig. 7. Layout only, data unchanged; Fig. 7 is unchanged pixel for pixel (`figures/fig6_order_comparison.py`, `figures/figstyle.py`).

## 9. Entry point for the nested validation

`code/nested_generic.py` holds the nested routines but has no command-line entry point, so the call
`python code/nested_generic.py` that the README listed among the reproduction commands ran without output. It is
replaced by `code/run_nested_observed.py`, which calls the functions the permutation scripts use for their observed
values (Route A1 `altorder.alt_nested_cv`, Route A1u `unsup_rep.unsup_nested_cv`, Route A2 `nested_generic.nested_cv`)
for all five strains and writes `results/recalc/nested_observed.json`. Run from a clean unpacking of this archive with
the Route A descriptor matrix in place, it reproduces all 30 honest Q2 and AUC values of Routes A1, A1u and A2 in
Table 1 (blocks a and b) to the two decimals printed there. No routine was changed.

## 10. Revised data files and stereocentre count in the reproduction path

`data/RouteB_descriptor_data.xlsx` is now the revised Additional file 3, byte for byte; the submitted file is kept
beside it as `data/RouteB_descriptor_data_submitted.xlsx`. Until now the bundled file was the submitted one, which the
README described as identical to the descriptor file of the article; that holds for the submitted article, not for the
revised one, whose Additional file 3 carries the corrected value. `code/desc_compute.py` counts stereocentres with the
CIP-based routine (`useLegacyImplementation=True`), as the revised file does (1,8-cineole 0 instead of 2), so
`code/add_structure_identifiers.py` verifies all 775 descriptor values of the revised file; with the submitted routine
it stopped at that value. On the revised file `run_all.py`, `fold_reduction.py`, `nlvcv.py` and `extra_tables.py`
reproduce the archived objects to within 3e-7 (the revised file stores values to six decimals), and bit for bit on the
submitted file; the tie-corrected p-values of section 7 are the only other differences. No reported number changes.

The same holds for the VIP workbook that `figures/figS2_figS3_vip.py` reads: `data/20260710_VIP_EffectSize_Tabellen.xlsx`
is now the revised Additional file 8, byte for byte, and the submitted workbook is kept as
`data/20260710_VIP_EffectSize_Tabellen_submitted.xlsx`. Until now the script redrew Figs. S3 and S4 as submitted from the
bundled workbook, and the revised versions were reachable only by pointing `QSAR_XLS` at the revised Additional file 8.
It now draws them as printed in the revised Supplement, pixel for pixel, and `figures/output/` holds these versions
(previously the submitted ones); with the submitted workbook it reproduces the submitted figures pixel for pixel.

Three gaps in the revision scripts are closed: `revision/make_base.py` rebuilds `base.pkl`, the
shared input of `revision/common.py` (verified identical to the file used, object for object); `revision/jobB_hits.py`
no longer reads a helper file of the working directory that is not part of this archive (its output is unchanged); and the instrumented module copies used by `revision/audit_fallback.py` and `revision/audit_gaucher.py` are
bundled as `revision/auditcode/` and `revision/auditcode_g/`.

## 11. Numbering of the additional files

For the revision the additional files are numbered in the order of their first citation in the manuscript. Every reference of the form "Additional file N" in this archive uses the new numbers, in the README files and in the comments and help texts of the scripts alike. File names keep the original numbers (for example `Additional_file_2_RouteB_descriptor_data.xlsx`), because the scripts open them by name, as do German comments with "Zusatzdatei"; only comments and documentation were rewritten, and for every script it was checked that its token stream without comments and strings is unchanged. Original number -> revised number: 1 -> 1, 2 -> 3, 3 -> 7, 4 -> 6, 5 -> 8, 6 -> 9, 7 -> 10, 8 -> 11, 9 -> 12, 10 -> 13; files new in the revision: 2, 4, 5.

## 12. Numbering of the supplementary figures and tables

For the revision the supplementary figures and tables are numbered in the order of their first citation in the main text. The text files of this archive use the new numbers. Scripts, figure file names (for example `figures/output/FigS2_VIP_per_descriptor.png`) and German comments written before the renumbering keep the original numbers so that the scripts still run. Original -> revised: Figs. S1 to S6 -> S2, S3, S4, S5, S6, S1; Tables S1 to S22 (S15 to S22 new in the revision) -> S1, S2, S6, S7, S8, S12, S14, S16, S17, S15, S18, S20, S21, S22, S4, S3, S5, S13, S9, S19, S10, S11.
