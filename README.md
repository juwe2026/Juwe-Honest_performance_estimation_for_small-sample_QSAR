# Honest performance estimation for small-sample QSAR

Analysis code and raw permutation outputs for:

> Werle J, Rondevaldova J, Werle E. *Honest performance estimation for small-sample QSAR: selection leakage, permutation-null inflation, and a label-free remedy.* Journal of Cheminformatics (submitted).

This repository is Additional file 3 of that article. It contains the analysis code and the
archived result objects (raw permutation outputs, cached intermediate results) that the numbers
reported in the manuscript and its Supplementary Information were computed from. See "Known
limitations" below before relying on `code/` to reproduce those numbers by re-running it.

Development version: https://github.com/juwe2026/Juwe-Honest_performance_estimation_for_small-sample_QSAR
Released under the MIT licence (see `LICENSE`). Each release is archived on Zenodo with its own
DOI; the article cites the DOI of the release that corresponds to the submitted version, and any
later release is reachable from the same Zenodo record.

## What the study does

Minimum inhibitory concentrations of 31 monoterpenoids against five pneumonia-associated
bacterial strains are dichotomised into five binary endpoints. Four descriptor-selection routes
are compared under apparent and under fully nested leave-one-out validation, with significance
assessed by permuting the class labels through the complete pipeline:

| Route | Reduction |
| --- | --- |
| A1 | correlation grouping first, group representative chosen by effect size, then supervised effect-size filter |
| A1u | as A1, but the representative is chosen by network centrality, which makes the grouping label-free |
| A2 | supervised effect-size filter first, correlation grouping of the survivors |
| B | 25 literature-derived descriptors reduced to 13 representatives by unsupervised correlation grouping; 12 of those are used for modelling, since one further pair (C_Count and EsterLacton_flag) is exactly redundant (r = rho = 1.000) and collapses to a single representative |

## Layout

```
code/                   analysis scripts (Python, RDKit, scikit-learn)
  manuscript_build/     scripts that assemble the manuscript and supplement documents
  add_structure_identifiers.py   writes the canonical SMILES, InChI and InChIKey columns
                         of Additional file 2 from its own SMILES column
data/                   Route B descriptor data; the Gaucher comparison data set
figures/                nine scripts generating Fig. 1-7 and Fig. S1-S6 (dataio.py and
                         figstyle.py are shared helpers, not figures themselves)
results/raw_permutations/   permutation outputs and per-route result objects
results/Additional_file_6_RouteA1_top15_EffSize_vs_VIP.xlsx
results/Additional_file_7_RouteA2_top15_EffSize_vs_VIP.xlsx
results/fold_reduction.json
results/fold_grouping_routeA.json
```

The Route A descriptor matrix (2,208 cleaned descriptors for 31 compounds) is not included
here because of its size; it is Additional file 4 of the article.

## Which file holds which result

| File | Contains |
| --- | --- |
| `final_all.npy` | Route A1, A1u, A2 and the top-15 variants: apparent and nested Q2/AUC with permutation p-values. Table 1, Table 3 and Table 5. |
| `vip1lv_final.npy` | Route A1 top-15 selected by VIP, one latent variable. Table 4. |
| `a2_top15_results.npy` | Route A2 top-15. Additional file 7. |
| `hits_three_variants.npy` | Observed and null univariate hit counts for A2, A1 and A1u, 1,000 permutations each. Table S6, Fig. 7, Fig. S6. |
| `routeB_final.json` | Route B with 13 representatives (superseded by the published 12-descriptor set; see `figures/data_recalc/routeB12.json`) and with the 25 a priori descriptors. Table 2 (25-descriptor rows only; the Table 1 Route B rows use the 12-descriptor set). |
| `routeB_nlv3.pkl`, `mono_nlv.pkl`, `mono_dcv.pkl` | Effect of the number of latent variables, including the doubly nested estimates. Table S5. |
| `lpo_all_routes.json`, `lpo_compare.pkl` | Leave-pair-out c-statistic and leave-2-out / leave-3-out Q2. Table S14. |
| `gaucher_res.pkl`, `gaucher_optlv.pkl`, `gaucher_doublecv.pkl`, `code/gaucher_doublecv.py` | The three routes applied to the SELDI-TOF data of Smit et al. (2007). Table S12; `gaucher_doublecv.py` regenerates the doubly nested row. |
| `gaucher_downsample.pkl`, `gaucher_bootstrap.pkl`, `code/gaucher_downsample.py`, `code/gaucher_bootstrap.py` | Detection-limit analysis by thinning the Gaucher data, and the bootstrap-the-controls check quoted in its legend. Table S13. The two `.pkl` files hold the earlier 150-draw runs. The table and its legend as published use the 2,000-draw runs, whose complete outputs are bundled as `results/raw_permutations/gaucher_downsample_2000.json` (all 10,001 per-draw Q2 and AUC values) and `results/raw_permutations/gaucher_bootstrap_2000.json` (both bootstrap arms, with the per-draw duplicate-control counts), so every figure in Table S13 and its legend can be recomputed from an archived object instead of repeating the run; `python code/gaucher_downsample.py 2000` and `python code/gaucher_bootstrap.py 2000` regenerate them. |
| `r2_vs_q2.json` | R2 against Q2 for Routes A1 and B. Table S10. |
| `ga_HI.pkl`, `ga_SA.pkl` | Genetic-algorithm wrapper comparison. Table S11. |
| `scaling_leakage.json` | Global versus in-fold autoscaling. |
| `fold_reduction.json` | Route B with the unsupervised correlation grouping recomputed inside every training fold, against the grouping fixed on all 31 compounds. Table S1. |
| `fold_grouping_routeA.json` | The same test for the 2,208-descriptor pool: Routes A1 and A1u with the 182 correlation groups rebuilt inside every training fold. Table S2. |
| `vip_simil.npy` | Spearman correlations between the per-strain VIP profiles, one latent variable. Superseded by the current 12-descriptor Route B set for some values (see `figures/data_recalc/routeB12.json`, `figures/data_recalc/vip_routeB.json`); not regenerated at 12 descriptors. |
| `s_two.json` | Additional two-stage VIP variant on the 25 literature descriptors. Table S8. |
| `perm*_*.npz` | The raw null distributions, 2,000 permutations per model, for Routes A1, A1u, A2 and the 13-/25-descriptor Route B variants. |
| `code/run_routeB12.py`, `results/raw_permutations/permB_12_*.npz`, `permB_12nlv3_*.npz` | The raw null distributions for the published 12-descriptor Route B model (one and three latent variables), 2,000 permutations per strain per model; reproduces `figures/data_recalc/routeB12.json`'s cached summary from scratch (observed Q2/AUC match to the last printed digit). Table 1, Table 2, Fig. 2, Fig. 4, Fig. 6. |
| `results/results.json` | Route B with both descriptor sets side by side, written by `code/run_all.py`: permutation tests at one and three components, VIP tables and confusion matrices. The `D12` block holds the published 12-descriptor values of Table 1 and Table 2; the `D13` block is the superseded 13-descriptor set, kept so the effect of dropping the redundant C_Count can be measured rather than asserted. |
| `figures/output/FigS3c_axis_key.csv` | The Fig. S3c x-axis position (1-60) to descriptor name key, written by `figures/figS2_figS3_vip.py`; republished as the `FigS3c_axis_key` sheet of Additional file 5. |

**Note.** The entry `A1v15` in `final_all.npy` is a superseded two-component variant of the
VIP top-15 model and does **not** correspond to Table 4. Use `vip1lv_final.npy` for that table.

The table references above match the manuscript's current Supplementary Information numbering
throughout.

## Reproducing

### One input has to be placed by hand

The Route A descriptor matrix is Additional file 4 of the article and is not bundled here
because of its size (17 MB). Download it and put it in `data/` under the name the code
expects:

```
data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx
```

In the submission package the same file is called
`Additional_file_4_RouteA_descriptor_matrix.xlsx`; the two are byte-identical, only the name
differs. To keep it elsewhere, point `QSAR_ROUTEA_XLSX` at it instead. Everything the Route B
and figure scripts need is already in `data/`.

Results are written to `results/recalc/`, or to `QSAR_RECALC_DIR` if that is set.

### The commands

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python code/desc_compute.py               # Route B descriptors from SMILES
python code/add_structure_identifiers.py data/RouteB_descriptor_data.xlsx
                                          # re-derives Additional file 2's canonical SMILES,
                                          #   InChI and InChIKey and checks them against the
                                          #   file; exits non-zero on any mismatch
python code/run_all.py                    # Route B, both descriptor sets: reproduces
                                          #   results/results.json (about 9 minutes)
python code/nested_generic.py             # nested validation, Routes A1 / A1u / A2
python code/fold_reduction.py 2000        # Table S1: grouping recomputed inside every fold
python code/fold_grouping_routeA.py 2000  # Table S2: the same for the 2,208 pool
python figures/fig1_scheme.py             # Fig. 1

# whole-pipeline permutations, one call per strain (Fig. 5, Table 3)
python code/prun_a1.py e HI    Cleaned_Steps1-4 2000 HI    100000
python code/prun_a1.py e SA    Cleaned_Steps1-4 2000 SA    100000
python code/prun_a1.py e SPneu Cleaned_Steps1-4 2000 SPneu 100000
python code/prun_a1.py e SPyo  Cleaned_Steps1-4 2000 SPyo  100000
python code/prun_a1.py e PA    Cleaned_Steps1-4 2000 PA    100000
```

`prun_a1.py` takes its arguments in the order `MODE SHORT SHEET NPERM CODE [BUDGET]`:

| Argument | Meaning |
| --- | --- |
| `MODE` | `e` selects the top 15 by effect size (Fig. 5, Table 3), `v` by VIP (Table 4) |
| `SHORT` | strain key: `HI`, `SA`, `SPneu`, `SPyo`, `PA`. It also selects the activity column |
| `SHEET` | sheet of Additional file 4 holding the descriptor pool: `Cleaned_Steps1-4` |
| `NPERM` | number of permutations; the published runs use 2,000 |
| `CODE` | file suffix of the output `.npz`, conventionally the same as `SHORT` |
| `BUDGET` | optional wall-clock limit in seconds. The run saves as it goes and resumes where it stopped, so a short budget can be called repeatedly |

`Cleaned_Steps1-4` is the pooled sheet and carries all five activity columns side by side, so
the strain has to be named rather than inferred. The mapping is:

| Strain | `SHORT` | Activity column in Additional file 4 |
| --- | --- | --- |
| *H. influenzae* | `HI` | `H. I. (H. influenzae)` |
| *S. aureus* | `SA` | `S. A. (S. aureus)` |
| *S. pneumoniae* | `SPneu` | `S. Pneu (S. pneumoniae)` |
| *S. pyogenes* | `SPyo` | `S. pyo (S. pyogenes)` |
| *P. aeruginosa* | `PA` | `P. A. (P. aeruginosa)` |

The same mapping is what `nested_generic.build(sheet, ycol=...)` expects. Its `ycol` defaults
to the third column, which is the convention of the strain-specific `CS1-5_*` sheets; on the
pooled sheet that default would run every strain against *H. influenzae*'s labels, so
`prun_a1.py` and `fold_grouping_routeA.py` set it explicitly.

### What it costs

Permutation runs are the expensive part: 2,000 permutations for five strains require roughly
0.3 million model fits under leave-one-out. In wall-clock terms that is less forbidding than it
sounds - about 130 seconds per strain, so eleven minutes for all five on one core. The
fold-wise regrouping of the 2,208-descriptor pool (Table S2) is the slower step at roughly five
seconds per fold.

### The Gaucher comparison data

```bash
python code/gaucher_downsample.py 2000    # Table S13: detection-limit analysis
python code/gaucher_bootstrap.py 2000     # Table S13 legend: the bootstrap check
python code/gaucher_doublecv.py           # Table S12, doubly nested row
python code/gperm.py A1 2000              # permutation null, one latent variable
python code/gperm2.py A1 2 2000           # permutation null, two latent variables
```

`gaucher_downsample.py` draws subsample `d` for case count `k` from
`numpy.random.default_rng([SEED, k, d])`, so every subsample is fixed by `(SEED, k, d)` alone
and any single draw can be regenerated on its own. `gaucher_doublecv.py` chooses the number of
latent variables by a second leave-one-out inside every training fold, over the candidates one
to five, and records the choice per fold; it takes about ten seconds per route.

## Conventions

- Every model uses one latent variable unless stated otherwise; the three-component Route B
  models are a sensitivity analysis, not the primary inference.
- `p = (1 + #{permuted >= observed}) / (1 + n_permutations)`.
- Centring and scaling parameters are estimated inside every fold, on the training compounds only.
- *P. aeruginosa* has three active compounds, is an artefact and is reported descriptively only.

`code/core.py` (the shared Route B pipeline) reads the bundled `data/RouteB_descriptor_data.xlsx`
via an overridable environment variable rather than a hardcoded, machine-specific path, uses the
current English sheet/column names, and uses nLV = 2 uniformly for the VIP table. The table
references and Route B description in this README match the manuscript's current numbering.
`figures/dataio.py` reads every path from an environment variable with a portable,
archive-relative default, and the same pattern is now used by every script on the documented
reproduction path - `run_all.py`, `desc_compute.py`, `nested_generic.py`, `altorder.py`,
`unsup_rep.py`, `prun_a1.py`, `fold_reduction.py`, `fold_grouping_routeA.py`,
`gaucher_downsample.py`, `gperm.py` and `gperm2.py` - each resolving inputs relative to the
archive root and writing to `results/recalc/` unless `QSAR_RECALC_DIR` overrides it.

`core.py` spells out both Route B descriptor sets explicitly. `RED12` is the 12 representatives
read from Additional file 2; `RED13` is those plus `C_Count`, the superseded set from before that
descriptor turned out to be perfectly collinear with `EsterLacton_flag`. Until now `RED13` was read
from the `VIP_reduced` sheet and `C_Count` filtered out of it, but that sheet has carried no
`C_Count` row since the reduction was documented in it, so the filter matched nothing and the two
sets were identical: `run_all.py` compared the 12-descriptor model with itself and would have
reported a difference of exactly zero. Re-running `run_all.py` now reproduces `results/results.json`
in full.

Two group caches, `altorder.get_groups()` and `unsup_rep.centrality_reps()`, key their stored
result on a hash of the descriptor matrix they were given. Earlier they wrote to a fixed file
name and returned it on any later call, so analysing a second sheet in the same working
directory silently reused the first sheet's groups; keying the cache makes a different matrix
use a different file.

`data/RouteB_descriptor_data.xlsx` is byte-identical to Additional file 2 of the article. It
was for a while an older generation of that file: the `Descriptors` sheet, which is the only one
any computation reads values from, was already identical, but its `VIP_reduced` sheet still
carried the superseded per-strain component numbers (4/4/3/1/4) instead of the uniform two that
Additional file 5, Fig. S2, Fig. S3 and `results/results.json` use, and its `Reduction` sheet
listed `EsterLacton_flag` twice, so `core.LIT25` came out with 26 entries rather than 25.
Neither affected a number: `core.py` reads only the descriptor names from `VIP_reduced`, and
`fold_reduction.py` de-duplicates `LIT25` before asserting its length. Both copies are now the
same file.

That file's `Descriptors` sheet carries three machine-readable structure identifiers next to the
input SMILES of each compound: the RDKit canonical SMILES, the standard InChI and the standard
InChIKey. `code/add_structure_identifiers.py` writes them, deriving all three from the SMILES
column itself, that is from the exact structure `code/desc_compute.py` computed the descriptors
from, and refusing to write unless it has first recomputed all 775 descriptor values from those
same SMILES and matched them against the sheet. The identifiers therefore describe the modelled
constitution: the structures carry no stereochemistry apart from the double-bond geometry given
for citral and geranic acid, while the compound names give the configuration of the material that
was tested. `core.py` addresses this sheet by column name, so the three added columns do not move
any value it reads; `run_all.py` reproduces `results/results.json` unchanged from the file with
them.

`data/` bundles the two source spreadsheets `figures/figS2_figS3_vip.py` and
`figures/figS4_figS5_scores.py` read via `dataio.xlsx()`
(`20260710_VIP_EffectSize_Tabellen.xlsx`, `20260802_FigS5_S6_Rohdaten.xlsx`), so Fig. S2, S3, S4
and S5 regenerate from this archive without any external file. All five figure-generation
scripts (`fig1_scheme.py`, `fig2_fig4_comparison.py`, `fig3_fig5_figS1_permutation.py`,
`fig6_order_comparison.py`, `figS4_figS5_scores.py`) pass `dpi=300` explicitly to `savefig()`
(via the shared `figstyle.save()` helper), so every PNG in `figures/output/` carries a pHYs
resolution chunk; without it most viewers assume 72 dpi. Fig. 6's strain labels are single-line at
a 30-degree rotation, its significance asterisks sit at one common height per bar group with a
glyph size chosen so adjacent marks stay visibly separated, and its axis label reads "descriptors
entering PLS-DA", matching the figure legend; Fig. 6's and Fig. 7's y-axes leave headroom above
the tallest bar; Fig. S5 is sized to fit the print area with a small margin. In Fig. 3, Fig. 5 and
Fig. S1, the bold strain-name title above each panel row is offset from the grey
permutations-below-axis note beneath it so the two no longer touch. `figures/figS2_figS3_vip.py`
sorts the descriptors shown in Fig. S2 and Fig. S3 by mean VIP with the descriptor's own name as
an explicit secondary key, so exact ties resolve the same way on every run; without it, Python's
per-process hash randomisation left the order of tied descriptors non-deterministic between runs,
which would have made Fig. S3c's axis numbering (and the position-to-descriptor key it also
writes to `figures/output/FigS3c_axis_key.csv`, republished as the `FigS3c_axis_key` sheet of
Additional file 5) unreproducible. The five compound names that differed between Additional
file 2/9 and Additional file 4 (`1,8-Cineole`, `3-Carene, α-Carene`, `m-Cymene`, `p-Cymene`,
`α-Pinene`) are now written the same way everywhere, including in the bundled
`data/RouteB_descriptor_data.xlsx` and `data/20260802_FigS5_S6_Rohdaten.xlsx`.

`results/raw_permutations/` holds `permB_12_*.npz` and `permB_12nlv3_*.npz`, the raw
2,000-permutation null distributions for the published 12-descriptor Route B model (one and
three latent variables), produced by `code/run_routeB12.py` and verified to reproduce
`figures/data_recalc/routeB12.json`'s cached observed Q2/AUC and p-values exactly, for every
route including the 12-descriptor Route B result.

The manuscript's availability statement is deliberately narrower than "every script here is a
reproduction path". It states that every reported result can be obtained either through the
documented path in "Reproducing" above, which has been run end to end from a clean unpack, or
from an archived raw result object in `results/`. It does not state that everything under
`code/` runs from this archive; "Known limitations" below names, file by file, what does not.

## Known limitations

The following are flagged here rather than silently left inconsistent:

- Every command in "Reproducing" above has been run end to end from a clean unpack of this
  archive with Additional file 4 placed as described, and reproduces the archived result objects
  exactly: the observed Q2 and AUC of all five strains, and the permutation draws element for
  element. What has *not* been verified end to end is the rest of `code/`: 30 of the 51 files
  there (26 Python and the 4 JavaScript files under `code/manuscript_build/`) still contain
  hardcoded absolute paths, mostly `/home/claude/...` output caches and `/mnt/user-data/...`
  input files. The list is not an estimate: it is `code/build_af2.py`, `build_excel.py`,
  `check_ccount_esterlacton.py`, `compute_vip.py`, `extra_tables.py`, `ga_run.py`,
  `lit_final.py`, `lit_nlv1.py`, `lit_plsda.py`, `nested_h.py`, `nested_write.py`, `nlvcv.py`,
  `perm_runner.py`, `prun.py`, `prun15.py`, `prun_a1_k.py`, `prun_a1u_vip.py`, `prun_alt.py`,
  `prun_uns.py`, `prun_vip.py`, `prun_vip1.py`, `recompute_a1_figs.py`, `reduction.py`,
  `run_ordinary_one.py`, `run_routeB.py`, `vip_fig.py` and the four files under
  `code/manuscript_build/`. None of them is called by the documented path above. These are working and repair scripts outside the documented path; every script
  the documented path actually calls has been made archive-relative. The paths
  originate from the working environment in which the code was developed; the use of a large
  language model in implementing this code is declared in the Methods section of the article,
  under "Use of large language models".
- The four scripts that first assembled the 25 literature descriptors (`build_excel.py`,
  `lit_plsda.py`, `lit_final.py`, `compute_vip.py`) do not run from this archive. They are the
  original working chain and expect intermediate objects that were never archived
  (`desc_table.pkl`, `desc_full.pkl`, `reduced_desc.npy`), they address the per-strain sheets of
  Additional file 4 under their former names `CS1-4_*` rather than the current `CS1-5_*`, and
  `lit_plsda.py` and `compute_vip.py` still use the former spelling `S. A. (Staph. Aureus)` of
  one activity column. The workbook `20260706_31MT_LitDeskriptoren_Analyse.xlsx` that these
  scripts write is an *output* of the chain, not an input it is missing. Its content is
  published: the descriptor values are Additional file 2 and the VIP and effect-size tables are
  Additional file 5. The same numbers are reachable and were re-derived through `core.py`, which
  does run from the bundled data, so no result of the paper depends on this chain; the gap is in
  the reproduction path of a superseded working script, not in the evidence.
- Two of the archived Gaucher result objects still have no generating script here:
  `gaucher_optlv.pkl` and `gaucher_top15.pkl`. Both are
  supporting analyses rather than sources of a published table row, and `gaucher_pipe.py`, the
  library they were built on, runs from the bundled data. `gaucher_res.pkl` is reconstructed
  exactly by that library; `code/gaucher_downsample.py` regenerates Table S13,
  `code/gaucher_bootstrap.py` the bootstrap check quoted in its legend and
  `code/gaucher_doublecv.py` the doubly nested row of Table S12, each verified against the
  archived object. The verification differs in strength and is worth stating exactly: the
  doubly nested script reproduces the archived Q2, AUC and the per-fold record of how many
  components were chosen to the last bit, for all three routes; the downsample script
  reproduces the deterministic k = 19 anchor to the last bit, and its sampled columns agree
  with the archived ones within the draw-to-draw variation, since the two runs use different
  random seeds; the bootstrap script likewise agrees distributionally and reproduces the mean
  number of duplicated controls. Where a table now reports the 2,000-draw figures, they come
  from these scripts and not from the archived 150-draw objects, which are kept for reference. The
  2,000-draw outputs themselves are bundled as
  `results/raw_permutations/gaucher_downsample_2000.json` and
  `results/raw_permutations/gaucher_bootstrap_2000.json`, so no figure in Table S13 or its legend
  depends on a run the reader has to repeat.
- Genuine German prose (docstrings, comments or console/log output, not just isolated column-name
  references) remains in a number of files under `code/`, concentrated in the working/repair
  scripts named below rather than the core pipeline (e.g. `s_two.json`'s own keys, "Erreger",
  "Q2_genestet", and full German sentences in `altorder.py`, `build_af2.py`,
  `fold_reduction.py`, `plsda.py`).
- `code/` mixes the analysis pipeline with working/repair scripts (`build_af2.py`,
  `fix_af2_af9.py`, `check_ccount_esterlacton.py`, `altorder.py`, `recompute_a1_figs.py`); it has
  not been reorganised into a clean, ordered pipeline. Two of them no longer run against the
  current data: `fix_af2_af9.py` expects a `Table1B_top15` sheet that the repair it performs has
  since renamed away, and `build_af2.py` rebuilds an earlier generation of Additional file 2.
  They are kept as a record of what was done, not as reproduction steps.

## Licence

MIT, see `LICENSE`.
