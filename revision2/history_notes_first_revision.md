# Notes on the first revision removed from the Additional files

For the second revision, the Additional files describe their content as it is and no longer how it changed. The notes
on the changes of the first revision (release 1.1.0) that were removed from them are kept here word for word, as they
stood in the Additional files of the first revision. Numbers inside the quotations follow the numbering of the first
revision. Each heading gives the published number of the Additional file and its number in the first revision.
`revision2/tables/aenderungen_af3_13.py` lists every change made to these files.

## Additional file 4 (first revision: Additional file 3), sheet Descriptors, AJ:AJ

The column `NumStereoCenters_submitted` (31 values, the stereocentre counts of the submitted file) was removed. The same
column is kept in `data/RouteB_descriptor_data.xlsx` of this archive.

## Additional file 4 (first revision: Additional file 3), sheet Computation_rules, E11

```
Assigned and unassigned CIP stereocentres; ring-bridgehead atoms that are not stereogenic are not counted (1,8-cineole = 0). The submitted version used useLegacyImplementation=False (1,8-cineole = 2); those values are kept in column NumStereoCenters_submitted of the Descriptors sheet.
```

## Additional file 4 (first revision: Additional file 3), sheet VIP_two_stage, A8

```
Revision: recomputed with the corrected stereocentre count of 1,8-cineole (0 instead of 2); algorithm and nLV unchanged.
```

## Additional file 4 (first revision: Additional file 3), sheet VIP_two_stage, A9

```
Submitted values (Q2 naive / nested): H.I. 0.287 / 0.161; S.A. 0.498 / 0.424; S.Pneu 0.525 / 0.305; S.Pyo 0.318 / 0.245; P.A. 0.998 / 0.995.
```

## Additional file 4 (first revision: Additional file 3), sheet VIP_two_stage, A10

```
The naive Q2 remains higher than the nested Q2 for every strain.
```

## Additional file 4 (first revision: Additional file 3), sheet Comparison_3_routes, A8

```
Regenerated for the revision from the result files behind Table 1 of the revised manuscript (Route A2 and Route B rows; apparent
```

## Additional file 4 (first revision: Additional file 3), sheet Comparison_3_routes, A9

```
values from the archived ordinary permutation runs). The sheet compares three conditions: the 2,208-descriptor pool of Route A2
```

## Additional file 4 (first revision: Additional file 3), sheet Comparison_3_routes, A12

```
p-values: 2,000 permutations, ties counted as in Table 1. The submitted version of this sheet contained values of an earlier stage.
```

## Additional file 4 (first revision: Additional file 3), sheet Conclusions, A1

```
CONCLUSIONS - literature-based descriptor selection (leakage-free); regenerated for the revision from the current results
```

## Additional file 4 (first revision: Additional file 3), sheet VIP_gt1_list, A8

```
Revision: recomputed with the corrected stereocentre count of 1,8-cineole (0 instead of 2).
```

## Additional file 4 (first revision: Additional file 3), sheet Info, A21:A39

```
REVISION (03.10.2026)
Software: RDKit 2026.03.6, Python 3.11.15. The RDKit version used for the submitted file was not recorded.
All 25 descriptor values were recomputed with this version from the input SMILES. Every submitted value was reproduced
except one: NumStereoCenters of 1,8-cineole, which the submitted file gave as 2. That file had used
FindMolChiralCenters(..., useLegacyImplementation=False), which also counts ring-bridgehead atoms as potential
stereo elements; the CIP-based count (useLegacyImplementation=True, the RDKit default) gives 0, the correct value for
this achiral molecule. The corrected value is used here and in the revised manuscript; the submitted values are kept
in column NumStereoCenters_submitted. Floating-point values are stored to six decimals.
Pearson_Corr and Spearman_Corr were recomputed from the stored values. Spearman coefficients use average ranks for ties.
In the submitted file they had been computed from unrounded values, in which isomers of the same formula differ in MolWt
and LabuteASA in the 14th decimal; this broke ties and caused the deviations of up to 0.019 noted in review.
Largest change of a Spearman coefficient against the submitted sheet: 0.2545 (NumStereoCenters / Ether_flag), including the effect of the
corrected stereocentre count. The 30 pairs with |r| or |rho| > 0.7 (25 by both coefficients, 3 by Spearman only, 2 by
Pearson only) and the reduction in the Reduction sheet are unchanged.
NumStereoCenters is not among the 12 modelled descriptors, so PLS-DA_reduced and VIP_reduced are unchanged apart from
p-values, which now count ties as in Table 1. VIP_two_stage and VIP_gt1_list, which use all 25 descriptors, were
recomputed with the corrected value.
The summary sheets Comparison_3_routes and Conclusions were regenerated from the current results (the submitted versions
contained values of an earlier stage); every number in them is read from the same result files as Table 1 and Table 2.
```

## Additional file 13 (first revision: Additional file 8), sheet README, A12

```
REVISION (03.10.2026): in the VIP_Lit sheets the columns TwoStage_VIP_apparent, TwoStage_VIP_nested_mean, TwoStage_selfreq_VIPgt1 and EffectSize, which use all 25 literature descriptors, were recomputed with the corrected stereocentre count of 1,8-cineole (0 instead of 2; see the Route B descriptor file). Rows are again sorted by TwoStage_VIP_apparent. The RouteB columns (12-descriptor set, which does not contain NumStereoCenters) and the ChemDes sheets are unchanged. The EffectSize of MolWt is now computed with ties between isomers of the same formula, which the submitted values had broken because of floating-point noise in the 14th decimal (change at most 0.024).
```

## Additional file 14 (first revision: Additional file 9), sheet Info, A24:A40

```
Revision of 06 October 2026: correction of the fallback prediction.
In the branch that applies when no descriptor passes the univariate filter within a fold, or when
the training fold contains only one class, the deposited code predicted the mean of the training
labels. Under leave-one-out that mean is (S - y_i)/(n - 1) and therefore carries the held-out label.
The fallback is now the neutral constant 0.5. The branch is never reached on the observed labels, in
any of the 31 folds, so every point estimate in this file is unchanged and reproduces exactly.
The permutation p-values have been recomputed with the corrected code, with exactly 2,000
permutations and with the seed of the published run, so that the difference is the effect of the
correction and not of a different permutation stream.
In the same recomputation, permuted values equal to the observed one are counted as at least as
extreme (tolerance 1e-9). AUC values are fractions k/(n1 n0) and equal fractions can be stored with
rounding differences of about 1e-16; the deposited code compared strictly and missed part of these
ties, which made some AUC p-values slightly too small.
Result: the fallback correction alone changes none of the sixteen p-values of the four evaluable
endpoints. Counting ties raises three AUC p-values of the effect-size arm (H. influenzae 0.0505 to
0.0555, S. pneumoniae 0.0195 to 0.0200, S. pyogenes 0.0035 to 0.0040); no verdict changes. Table 3
and Table 4 of the main text carry the same values.
```

## Additional file 15 (first revision: Additional file 10), sheet Info, A14

```
Permutations of the whole pipeline: about 2,000 per model (strict parity).
```

## Additional file 15 (first revision: Additional file 10), sheet Info, A22

```
only for S. pneumoniae and S. pyogenes. For S. aureus its AUC p-value is 0.0525 (previously 0.0497,
```

## Additional file 15 (first revision: Additional file 10), sheet Info, A23

```
see the revision note below) and its Q2 p-value 0.060, so neither reaches 0.05 under effect size.
```

## Additional file 15 (first revision: Additional file 10), sheet Info, A30:A59

```
Revision of 06 October 2026: correction of the fallback prediction.
In the branch that applies when no descriptor passes the univariate filter within a fold, or when
the training fold contains only one class, the deposited code predicted the mean of the training
labels. Under leave-one-out that mean is (S - y_i)/(n - 1) and therefore carries the held-out label.
The fallback is now the neutral constant 0.5. The branch is never reached on the observed labels, in
any of the 31 folds, so every point estimate in this file is unchanged and reproduces exactly.
The permutation p-values have been recomputed with the corrected code, with exactly 2,000
permutations and with the seed of the published run, so that the difference is the effect of the
correction and not of a different permutation stream.
In the same recomputation, permuted values equal to the observed one are counted as at least as
extreme (tolerance 1e-9). AUC values are fractions k/(n1 n0) and equal fractions can be stored with
rounding differences of about 1e-16; the deposited code compared strictly and missed part of these
ties, which made some AUC p-values slightly too small.
Result: the fallback correction alone changes p-values by at most 0.0011. Counting ties raises
several AUC p-values by up to 0.003; one nominal verdict changes: the effect-size arm for S. aureus
has an AUC p-value of 0.0525 instead of 0.0497 and is no longer below 0.05 (it lies within Monte
Carlo error of 0.05). The deposited runs of this file had
1,900 to 2,000 permutations rather than the nominal 2,000, because permutations with a degenerate
training fold were dropped. Part of the small differences is that change in count.

Correction to this Info sheet: the VIP ranking of the Route A2 arm is computed from a
two-component PLS-DA, not a one-component model as previously stated here. The archived values
reproduce exactly with two components and differ by up to 0.11 in Q2 with one. The cause is
visible in the code: the call omits the nlv_vip argument and the default is two. The Route A1 arm
in Additional file 9 does use one component, as stated there.

The function this arm was computed with, nested_cv_viptop, was missing from the deposited
archive. It has been reconstructed and reproduces all ten archived values exactly. The script
code/prun_vip.py now calls the reconstruction and names the activity column of each strain
explicitly (CHANGELOG_revision.md, section 5).
```

## Additional file 16 (first revision: Additional file 12), sheet Info, A9:A13

```
LV1_score in FigS5_LV1_scores. For H. influenzae it did not (Component1 was exactly the negative of
LV1_score, e.g. (-)-Menthone: -2.66055 in the original computation vs. +2.66055 for LV1_score); both
Component1 in FigS6_scores_2comp and Loading_Comp1 in S6_loadings_H.I. have been sign-flipped here so
that all five strains follow the same convention as Fig. S5 (active compounds score higher on
Component1). Component2 and the autoscaled descriptor matrices are unaffected by this change.
```

## Additional file 7 (first revision: Additional file 13), sheet Info, A14:A16

```
Revision (October 2026): the p-values in sheet Routes_A1_A1u_A2 are now the exact permutation p-values from the archived
null distributions (2,000 permutations; p = (1 + number of permuted values >= observed) / 2,001). 0.0005 is the smallest attainable
value and means that no permuted value reached the observed one. The earlier entry "<0.001" was not exact for every route.
```
