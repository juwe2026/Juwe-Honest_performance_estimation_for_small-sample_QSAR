# The four parts of this archive in the journal submission

The journal accepts additional files of at most 20 MB. The complete archive (release 1.2.0, Zenodo) is therefore
submitted in four parts. Each part is a zip file with the same top folder. Unpacking all four into one folder gives the
complete archive.

| Part | Additional file | Content |
| --- | --- | --- |
| 1 | 10 | Code of all analyses, the data and results of the data sets of this study, the figures, the benchmark scripts, specifications, workbooks, control-compound tables and evaluated results (`benchmark/ergebnisse*/*.json`), and this file |
| 2 | 19 | Main benchmark analysis: `benchmark/datasets_k/` (derived data sets) and `benchmark/ergebnisse_k/perm/` (raw permutation outputs) |
| 3 | 20 | Benchmark analysis under the original rule: `benchmark/datasets/` and `benchmark/ergebnisse/perm/` |
| 4 | 21 | Median-rule sensitivity analysis: `benchmark/datasets_m/` and `benchmark/ergebnisse_m/perm/`, and the candidate lists `benchmark/cand_all.json` and `benchmark/cand_ranked.json` |

Part 1 alone reproduces every result of the data sets of this study. The benchmark evaluation also needs the data sets
and raw permutation outputs of parts 2 to 4. The benchmark data sets in parts 2 to 4 are
distributed under CC BY-SA 3.0 (`benchmark/DATA_LICENCE.md`).
