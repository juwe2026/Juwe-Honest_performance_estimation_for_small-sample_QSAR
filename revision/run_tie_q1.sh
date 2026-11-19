#!/bin/bash
# Archive-relative: runs in revision/ of this archive. QSAR_ROUTEA_XLSX defaults to the place named in README.md.
cd "$(dirname "$0")"
ROOT="$(cd .. && pwd)"
export QSAR_ROUTEA_XLSX="${QSAR_ROUTEA_XLSX:-$ROOT/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx}"
export RESUME=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python jobI2_fallback.py > logT_I2.txt 2>&1
python jobL_a2hi.py > logT_L.txt 2>&1
python jobJ2_topk.py > logT_J2.txt 2>&1
cd lone && python run3.py C > logT_run3C.txt 2>&1
cd "$ROOT/revision" && echo Q1 fertig > logT_q1_done.txt
