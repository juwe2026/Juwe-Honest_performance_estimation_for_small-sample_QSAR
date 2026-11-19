# -*- coding: utf-8 -*-
"""Proposed changes to the submitted Additional files 3, 4, 6 and 8 to 13 (P3, option A, EW 06.11.2026 13:56 Berlin).

Old number -> new number (approved N2): 3 -> 4, 4 -> 6, 6 -> 3, 8 -> 13, 9 -> 14, 10 -> 15, 11 -> 11, 12 -> 16, 13 -> 7.
Additional file 5 stays unchanged (document of record).

Kinds of change: T title and numbers, H reference to an earlier version removed or rewritten in the present tense,
S semicolon replaced, K content brought in line with the main text (corrected threshold 0.0167 for three effective
tests, Table S12), M mechanical renumbering (computed by apply_af3_13.py, not listed here).
Texts given here are final, already in the new numbering. Every other text cell is renumbered mechanically
(nummern_freigegeben.rn_old_text)."""

FILES = {
    3: dict(new=4, src='Additional_file_3_RouteB_descriptor_data 06.10.2026.xlsx',
            name='Additional_file_4_RouteB_descriptor_data'),
    4: dict(new=6, src='Additional_file_4_second_MIC_data_set 06.10.2026.xlsx',
            name='Additional_file_6_second_MIC_data_set'),
    6: dict(new=3, src='Additional_file_6_RouteA_descriptor_matrix 06.10.2026.xlsx',
            name='Additional_file_3_RouteA_descriptor_matrix'),
    8: dict(new=13, src='Additional_file_8_VIP_and_effect_size_tables 06.10.2026.xlsx',
            name='Additional_file_13_VIP_and_effect_size_tables'),
    9: dict(new=14, src='Additional_file_9_RouteA1_top15_EffSize_vs_VIP 06.10.2026.xlsx',
            name='Additional_file_14_RouteA1_top15_EffSize_vs_VIP'),
    10: dict(new=15, src='Additional_file_10_RouteA2_top15_EffSize_vs_VIP 06.10.2026.xlsx',
             name='Additional_file_15_RouteA2_top15_EffSize_vs_VIP'),
    11: dict(new=11, src='Additional_file_11_classification_metrics 06.10.2026.xlsx',
             name='Additional_file_11_classification_metrics'),
    12: dict(new=16, src='Additional_file_12_source_data_supplementary_figures 06.10.2026.xlsx',
             name='Additional_file_16_source_data_supplementary_figures'),
    13: dict(new=7, src='Additional_file_13_Gaucher_reanalysis 06.10.2026.xlsx',
             name='Additional_file_7_Gaucher_reanalysis'),
}

# (old file, sheet, cell) -> (kind, new text). None deletes the content of the cell.
CELLS = {
    # ---------------------------------------------------------------- old 3 -> new 4
    (3, 'Info', 'A15'): ('S', 'double-bond geometry specified for Citral and Geranic acid. The compound names give the configuration of'),
    (3, 'Computation_rules', 'E11'): ('H S', 'Assigned and unassigned CIP stereocentres. Ring-bridgehead atoms that are not stereogenic are not counted (1,8-cineole = 0).'),
    (3, 'Computation_rules', 'E20'): ('S', 'Ether, ester oxygen excluded'),
    (3, 'Computation_rules', 'E24'): ('S', 'Ertl method, N and O only'),
    (3, 'Reduction', 'A23'): ('S', 'across the remaining C10 series. The modelling set therefore comprises 12 descriptors, not the 13 shown'),
    (3, 'VIP_reduced', 'A15'): ('S M', "nLV = 2 for every strain (apparent VIP, matching Fig. S2b/Fig. S3a and Additional file 13's convention). The footnote below gives its own source rather than pointing to VIP_gt1_list (25-descriptor set). C_Count is not listed: it is exactly redundant with EsterLacton_flag (r = rho = 1.000, see the Reduction sheet) and is excluded from the 12-descriptor modelling set, so it has no VIP in this set. Source: figures/data_recalc/vip_routeB.json, key D12, field 'apparent'."),
    (3, 'VIP_two_stage', 'A8'): ('H', None),
    (3, 'VIP_two_stage', 'A9'): ('H', None),
    (3, 'VIP_two_stage', 'A10'): ('H', 'The naive Q2 is higher than the nested Q2 for every strain.'),
    (3, 'Comparison_3_routes', 'G2'): ('K S', 'Route B: significant after correction for three effective tests (alpha = 0.0167). Route A2 nested: not significant'),
    (3, 'Comparison_3_routes', 'G3'): ('K S', 'Route B: significant after correction for three effective tests (alpha = 0.0167). Route A2 nested: significant after the same correction (p = 0.0130, within Monte Carlo error of the threshold)'),
    (3, 'Comparison_3_routes', 'G4'): ('K S', 'Route B: significant after correction for three effective tests (alpha = 0.0167). Route A2 nested: significant after the same correction'),
    (3, 'Comparison_3_routes', 'G5'): ('S', 'Route B: nominally significant (p < 0.05). Route A2 nested: not significant'),
    (3, 'Comparison_3_routes', 'A8'): ('H S', 'Values from the result files behind Table 1 (Route A2 and Route B rows), apparent'),
    (3, 'Comparison_3_routes', 'A9'): ('H', 'values from the archived ordinary permutation runs. The sheet compares three conditions: the 2,208-descriptor pool of Route A2'),
    (3, 'Comparison_3_routes', 'A12'): ('H', 'p-values: 2,000 permutations, ties counted as in Table 1.'),
    (3, 'Conclusions', 'A1'): ('H S', 'CONCLUSIONS - literature-based descriptor selection (leakage-free)'),
    (3, 'Conclusions', 'A5'): ('K', '   Q2 = 0.29, 0.39 and 0.35, p = 0.0020, 0.0010 and 0.0015, all below the corrected threshold of 0.0167 for three effective tests (Table 1).'),
    (3, 'Conclusions', 'A7'): ('S', '   0.15 for S. aureus (p = 0.0130) and 0.19 for S. pneumoniae (p = 0.0095), see sheet Comparison_3_routes.'),
    (3, 'Conclusions', 'A9'): ('S', '   supervised selection from many candidates at n = 31. A small, theory-based, leakage-free set makes the moderate signal visible.'),
    (3, 'Conclusions', 'A12'): ('S', '4) S. pyogenes remains weak with the 12 descriptors (Q2 = 0.05, p = 0.0360, nominal only, AUC p = 0.1219, not significant), despite the'),
    (3, 'Conclusions', 'A13'): ('S M', '   largest active class (11). With all 25 descriptors it reaches Q2 = 0.24 (p = 0.0040, Table S16): the reduction removed descriptors'),
    (3, 'Conclusions', 'A16'): ('K', '6) Multiplicity: against the corrected threshold of 0.0167 for three effective tests (Table S12), the Q2 of H. influenzae, S. aureus and S. pneumoniae'),
    (3, 'Conclusions', 'A17'): ('K S', '   passes. On AUC, H. influenzae (p = 0.0015) and S. pneumoniae (p = 0.0075) pass and S. aureus (p = 0.0175) does not, within Monte Carlo error of the threshold.'),
    (3, 'Conclusions', 'A21'): ('S', '  cross-validation fold (nested). Otherwise it is the same leakage as with the univariate filter, with VIP as the criterion.'),
    (3, 'Conclusions', 'A27'): ('S', 'nLV that maximises Q2 on all compounds contain the optimism of that choice and are an upper bound. The selection-adjusted'),
    (3, 'Conclusions', 'A33'): ('S', '- Pearson alone gives 14 representatives instead of 13 (before resolving the exactly redundant pair, afterwards 13 and 12):'),
    (3, 'VIP_gt1_list', 'A8'): ('H', None),
    # ---------------------------------------------------------------- old 4 -> new 6
    (4, 'Info', 'B3'): ('S M', 'Lone A, Thomsen TT, Nielsen JE, Thulstrup PW, Klitgaard RN, Lobner-Olesen A, Lund R, Jenssen H, Hansen PR (2019) Structure-activity study of an all-D antimicrobial octapeptide D2D. Molecules 24:4571 (reference 25 of the main text). Published under CC BY 4.0. MIC values and sequences are reproduced from Table 1 of that article.'),
    (4, 'Info', 'B4'): ('S', 'Broth microdilution in Mueller-Hinton II broth in triplicate, two-fold dilutions from 128 to 2 ug/mL. ">128" means no inhibition at the highest concentration tested and is right-censored.'),
    (4, 'Info', 'B5'): ('S', 'SMILES reconstructed from the published sequences with RDKit. All residues are in D configuration (checked: every chiral alpha-carbon is R, plus the S beta-carbon of D-threonine in peptides 14 and 22). Calculated molecular weights agree with the published molecular weights (Supplementary Table S1 of the source) within 0.01 g/mol for all 37 peptides.'),
    (4, 'Info', 'B6'): ('S', 'Main analysis (rule A): active if MIC <= 8 ug/mL, for P. aeruginosa <= 16 ug/mL. Sensitivity analysis (rule C): active if MIC <= 16 ug/mL for all strains. Both rules were fixed in writing before any model was fitted (Additional file 5). Class sizes, rule A: E. coli 20, S. aureus 23, P. aeruginosa 14, A. baumannii 14 active of 37, rule C: 24, 26, 14, 24.'),
    (4, 'Info', 'B7'): ('S', '1,830 two-dimensional descriptors (Mordred 1,613, RDKit 217) were computed from the SMILES. After removal of non-numeric columns, columns with more than 50 % missing values and constant columns, and median imputation of the remaining gaps, 1,275 descriptors remain (sheet Descriptors_1275).'),
    # ---------------------------------------------------------------- old 6 -> new 3
    (6, 'Info', 'B14'): ('S', 'Columns with >50% missing values removed, remaining gaps median-imputed'),
    (6, 'Info', 'C14'): ('S', 'threshold 50%, 362 columns imputed'),
    (6, 'Info', 'A31'): ('S', "'1,8-Cineole' (not '1,8-cineol'), '3-Carene, α-Carene' (not '3-Carene. α-Carene', where a stray full"),
    (6, 'Info', 'A32'): ('S', "stop made two related isomers read as one entry), and 'm-Cymene', 'p-Cymene', 'α-Pinene' with the"),
    (6, 'Info', 'A35'): ('S M', "((-)-Menthone, (+)-Dihydrocarvone, (S)-(-)-Perillaldehyde, alpha,beta-Thujone), consistent with Table S7's"),
    (6, 'Apparent_PLSDA_results', 'A8'): ('S M', 'Apparent = one PLS-DA with one latent variable fitted once on all 31 compounds using the FIN_Des_* descriptor set of the corresponding strain (Route A2: univariate filter first, then correlation grouping of the survivors). Not cross-validated, reported here because the Additional file 3 description promises apparent PLS-DA results.'),
    # ---------------------------------------------------------------- old 8 -> new 13
    (8, 'README', 'A8'): ('S M', 'Route B uses 12 of the 25 literature descriptors (unsupervised correlation reduction, one further pair is exactly redundant and reduces to a single representative, see Table S8). Empty RouteB cells = removed in the correlation step.'),
    (8, 'README', 'A9'): ('S', 'Note: VIP here is from two components (as in MetaboAnalyst). The validated Q² models used one component.'),
    (8, 'README', 'A10'): ('S M', 'FigS3c_axis_key: the Fig. S3c x-axis is numbered 1-60 (union of the five strain-specific Route A1 top-15 sets, ordered by mean VIP). This sheet gives the descriptor name for each position number.'),
    (8, 'README', 'A12'): ('H S', 'In the VIP_Lit sheets the columns TwoStage_VIP_apparent, TwoStage_VIP_nested_mean, TwoStage_selfreq_VIPgt1 and EffectSize use all 25 literature descriptors, with the CIP-based stereocentre count of 1,8-cineole (0, see Additional file 4). Rows are sorted by TwoStage_VIP_apparent. The RouteB columns use the 12-descriptor set, which does not contain NumStereoCenters. The EffectSize of MolWt treats isomers of the same formula as ties.'),
    # ---------------------------------------------------------------- old 9 -> new 14
    (9, 'Info', 'A5'): ('S', 'the group member with the largest effect size. The supervised univariate filter is then'),
    (9, 'Info', 'A10'): ('S', 'on the Route A1 set. A second one-component PLS-DA on those 15 is the reported model.'),
    (9, 'Info', 'A13'): ('S', 'app = selection performed once on all 31 compounds, nested = selection repeated inside every fold.'),
    (9, 'Info', 'A16'): ('T', 'These are the underlying values of Table S17 (EffSize) and Table S18 (VIP) in Additional file 1.'),
    (9, 'Info', 'A19'): ('S', 'P. aeruginosa has three active compounds, is an artefact and is reported descriptively only.'),
    (9, 'Info', 'A20'): ('S', 'Its permutation p-values are therefore not given.'),
    # ---------------------------------------------------------------- old 10 -> new 15
    (10, 'Info', 'A3'): ('S', 'ROUTE A2. The supervised univariate filter is applied first, to all 2,208 cleaned descriptors.'),
    (10, 'Info', 'A4'): ('S', 'The survivors are then grouped by unsupervised correlation. The 15 highest-ranked descriptors'),
    (10, 'Info', 'A7'): ('T', 'These values are NOT the values of Tables S17 and S18 in Additional file 1. Those tables report Route A1,'),
    (10, 'Info', 'A8'): ('T', 'in which the grouping precedes the filter, and their underlying values are in Additional file 14.'),
    (10, 'Info', 'A13'): ('S', 'app = selection performed once on all 31 compounds, nested = selection repeated inside every fold.'),
    (10, 'Info', 'A14'): ('H', 'Permutations of the whole pipeline: 2,000 per model.'),
    (10, 'Info', 'A19'): ('S', '(effect size: S. pneumoniae and S. pyogenes, VIP: S. aureus and S. pyogenes). S. aureus and'),
    (10, 'Info', 'A20'): ('S', 'S. pneumoniae exchange significance for Q2. Only S. pyogenes is robustly significant under both.'),
    (10, 'Info', 'A22'): ('H', 'only for S. pneumoniae and S. pyogenes. For S. aureus its AUC p-value is 0.0525, within Monte Carlo error'),
    (10, 'Info', 'A23'): ('H', 'of 0.05, and its Q2 p-value 0.060, so neither reaches 0.05 under effect size.'),
    (10, 'Info', 'A25'): ('S', 'Under Route A1, at the same model size, the two criteria agree on every significance verdict.'),
    (10, 'Info', 'A26'): ('K', 'This difference between the two orders is discussed in Additional file 1 (section on fixed model size and VIP-based selection).'),
    # ---------------------------------------------------------------- old 11 -> new 11
    (11, 'Info', 'A11'): ('S M', 'representative, see Table S8). The confusion-matrix counts themselves are unchanged by that'),
    (11, 'Info', 'A13'): ('S', '13-descriptor sets). This equivalence has been checked only for leave-one-out cross-validation'),
    # ---------------------------------------------------------------- old 13 -> new 7
    (13, 'Info', 'A3'): ('M', 'to the SELDI-TOF proteomics data of Smit and colleagues (2007) [27], together with the 15 highest-'),
    (13, 'Methods', 'A3'): ('S', 'Input: Gaucherdata_paper39_median_normalised.csv. Autoscaling is performed WITHIN every cross-validation fold.'),
    (13, 'Methods', 'A8'): ('S', 'apparent = selection performed once on all samples, nested = selection repeated fully inside every fold.'),
}

# (old file, sheet): (first row, last row, kind, new lines). The rows first..last are replaced by the new lines,
# surplus rows are deleted. Column A, text only.
BLOCKS = {
    (3, 'Info'): (21, 39, 'H S', [
        'Software: RDKit 2026.03.6, Python 3.11.15. All 25 descriptor values were computed with this version from the input SMILES.',
        'NumStereoCenters is the CIP-based count (FindMolChiralCenters with useLegacyImplementation=True, the RDKit default). It does',
        'not count ring-bridgehead atoms that are not stereogenic and gives 0 for the achiral 1,8-cineole.',
        'Floating-point values are stored to six decimals.',
        'Pearson_Corr and Spearman_Corr are computed from the stored values, so that isomers of the same formula, which have the',
        'same MolWt and LabuteASA, are ties. Spearman coefficients use average ranks for ties.',
        'The 30 pairs with |r| or |rho| > 0.7 (25 by both coefficients, 3 by Spearman only, 2 by Pearson only) define the',
        'reduction in the Reduction sheet. NumStereoCenters is not among the 12 modelled descriptors. VIP_two_stage and',
        'VIP_gt1_list use all 25 descriptors. The p-values count ties as in Table 1.',
        'Every number in the summary sheets Comparison_3_routes and Conclusions is read from the same result files as',
        'Table 1 and Table S16 (Additional file 1).']),
    (9, 'Info'): (24, 40, 'H', [
        'Within a fold in which no descriptor passes the univariate filter, or in which the training fold contains',
        'only one class, the model predicts the neutral constant 0.5. This branch is never reached on the observed labels.',
        'Permuted values that equal the observed one within a tolerance of 1e-9 are counted as at least as extreme.']),
    (10, 'Info'): (30, 59, 'H', [
        'Within a fold in which no descriptor passes the univariate filter, or in which the training fold contains',
        'only one class, the model predicts the neutral constant 0.5. This branch is never reached on the observed labels.',
        'Permuted values that equal the observed one within a tolerance of 1e-9 are counted as at least as extreme.',
        '',
        'The VIP ranking of the Route A2 arm is computed from a two-component PLS-DA (the call does not set the nlv_vip',
        'argument, whose default is two). The Route A1 arm in Additional file 14 uses one component.',
        'The function nested_cv_viptop that computes this arm is in code/a2vip.py of the code archive (Additional file 10)',
        'and reproduces all ten values of this file exactly.']),
    (12, 'Info'): (9, 13, 'H S M', [
        'LV1_score in FigS7_LV1_scores. For H. influenzae the fit returned Component1 with the opposite sign (for example',
        '(-)-Menthone -2.66055 against +2.66055 for LV1_score). Component1 in FigS8_scores_2comp and Loading_Comp1 in',
        'S8_loadings_H.I. are therefore sign-flipped, so that all five strains follow the same convention as Fig. S7',
        '(active compounds score higher on Component1). Component2 and the autoscaled descriptor matrices are unaffected.']),
    (13, 'Info'): (14, 16, 'H', [
        'The p-values in sheet Routes_A1_A1u_A2 are the exact permutation p-values from the archived null distributions',
        '(2,000 permutations, p = (1 + number of permuted values >= observed) / 2,001). 0.0005 is the smallest attainable',
        'value and means that no permuted value reached the observed one.']),
}

# column with the values of an earlier version, deleted (old file, sheet, column letter, header)
DELETE_COLUMNS = [(3, 'Descriptors', 'AJ', 'NumStereoCenters_submitted')]

# sheet names that carry figure numbers, and their mentions in text cells
SHEET_RENAMES = {
    8: {'FigS4c_axis_key': 'FigS3c_axis_key'},
    12: {'FigS5_LV1_scores': 'FigS7_LV1_scores', 'FigS6_scores_2comp': 'FigS8_scores_2comp',
         'S6_loadings_': 'S8_loadings_', 'S6_matrix_': 'S8_matrix_'},          # prefixes of five sheets each
}
# the same sheet name in Additional file 1 (Fig. S3 legend)
AF1_SHEETNAME = ('(sheet FigS4c_axis_key)', '(sheet FigS3c_axis_key)')

# semicolons in the five sheets Info_Step5_<strain> of old file 6: '; ' -> ', ' in these cells
SEMI_TO_COMMA = [(6, 'Info_Step5_', ['A4', 'A6', 'A7', 'A10'])]
