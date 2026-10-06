#!/bin/bash
cd /home/claude/rev; export QSAR_ROUTEA_XLSX=/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx; export RESUME=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python jobI2_fallback.py > logT_I2.txt 2>&1
python jobL_a2hi.py > logT_L.txt 2>&1
python jobJ2_topk.py > logT_J2.txt 2>&1
cd lone && python run3.py C > logT_run3C.txt 2>&1
echo Q1 fertig > /home/claude/rev/logT_q1_done.txt
