#!/bin/bash
cd /home/claude/rev; export QSAR_ROUTEA_XLSX=/home/claude/rev/QSAR_selection_leakage/data/20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx; export RESUME=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python jobA_nofilter.py > logT_A.txt 2>&1
python jobC_routeB.py > logT_C.txt 2>&1
python jobE_misc.py > logT_E.txt 2>&1
python jobE2_morgan_fixed.py > logT_E2.txt 2>&1
python jobD_rf.py > logT_D.txt 2>&1
python jobD2_rf_morgan.py > logT_D2.txt 2>&1
python jobH2_median_fixed.py HI > logT_H2HI.txt 2>&1
python jobH2_median_fixed.py SPyo > logT_H2SPyo.txt 2>&1
python jobF2_no3d_fixed.py > logT_F2.txt 2>&1
cd lone && python run2_tie.py A > logT_run2A.txt 2>&1 && python run2_tie.py C > logT_run2C.txt 2>&1
python run3.py A > logT_run3A.txt 2>&1
echo Q2 fertig > /home/claude/rev/logT_q2_done.txt
