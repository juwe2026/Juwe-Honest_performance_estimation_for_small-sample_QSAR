# Licence of the ChEMBL-derived data in this folder

The ChEMBL benchmark uses data from ChEMBL release 37 (database version of 1 May 2026, released on 29 May 2026),
European Bioinformatics Institute (EMBL-EBI), https://www.ebi.ac.uk/chembl/. ChEMBL data are made available under the
Creative Commons Attribution-ShareAlike 3.0 Unported licence (CC BY-SA 3.0,
https://creativecommons.org/licenses/by-sa/3.0/).

Reference: Zdrazil B, Felix E, Hunter F, Manners EJ, Blackshaw J, Corbett S et al (2024) The ChEMBL Database in 2023:
a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Res 52(D1):D1180-D1192.
https://doi.org/10.1093/nar/gkad1004

The data files in this folder that contain ChEMBL records or values derived from them are distributed under the same
licence, CC BY-SA 3.0. These are the candidate and selection files (`cand_all.json`, `cand_ranked.json`,
`auswahl14.json`, `auswahl16.json`, `elternpools.json`, `doi2chembl16.json`, `free.csv`), the data sets in
`datasets/`, `datasets_k/` and `datasets_m/` (structures, ChEMBL identifiers, activity values, classes and the descriptor
matrices computed from the structures), the files in `kontrollen/`, the result files in `ergebnisse/`, `ergebnisse_k/`
and `ergebnisse_m/` and the results workbooks in `workbooks/`. The changes made to the ChEMBL records are described in `README_benchmark.md`
and in the specifications in `specifications/`.

The database file the selection scripts query (`chembl37_subset.db`, 2.1 GB) is not included. It holds the tables and
columns listed in `stream_parse.py`, which builds it from the ChEMBL 37 MySQL dump read on standard input.

The code in this folder is covered by the MIT licence of the archive (`../LICENSE`).
