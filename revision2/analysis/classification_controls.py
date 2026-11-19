# -*- coding: utf-8 -*-
"""Classification metrics of the control arms of Table 1 and Table S14 (response to Reviewer 2, point 6).

Decision threshold 0.5, a prediction of 0.5 or more counts as active (as in revision/common.py, mcc_conf).
  PLS-DA on all 2,208 descriptors and on the 182 label-free representatives: archived leave-one-out predictions
      (results/revision/bootstrap_predictions_ds1.json), checked against the counts in A_nofilter.json
  random forests (all 2,208, 182 representatives, Route B 12, Morgan bits): out-of-bag predictions, refitted
      exactly as in revision/jobD_rf.py (300 trees, square-root feature sampling, random_state 0); the
      out-of-bag Q2 and AUC must equal those in D_rf.json and D_rf_morgan.json, and for all 2,208 the
      predictions must equal the archived ones
  PLS-DA on Morgan bits, without selection and with the nested filter: archived counts (morgan2_fixed.json)
Usage: python classification_controls.py   (writes classification_controls.json to QSAR_REV2_OUT, default revision2/results)
"""
import os as _osp, sys
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_OUT = _osp.environ.get('QSAR_REV2_OUT', _osp.path.join(_ROOT, 'revision2', 'results'))
sys.dont_write_bytecode = True
sys.path.insert(0, _osp.path.join(_ROOT, 'code'))
import json
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import core, nested_generic as ng, altorder as ao, unsup_rep as ur

R = _osp.path.join(_ROOT, 'results', 'revision')
P = json.load(open(_osp.path.join(R, 'bootstrap_predictions_ds1.json')))
NF = json.load(open(_osp.path.join(R, 'A_nofilter.json')))
DRF = json.load(open(_osp.path.join(R, 'D_rf.json')))
MF = json.load(open(_osp.path.join(R, 'morgan2_fixed.json')))
MB = np.load(_osp.path.join(R, 'morgan2_bits.npy'))
KEYS = ['HI', 'SA', 'SPneu', 'SPyo', 'PA']
YCOL = {'HI': 'H. I. (H. influenzae)', 'SA': 'S. A. (S. aureus)', 'SPneu': 'S. Pneu (S. pneumoniae)',
        'SPyo': 'S. pyo (S. pyogenes)', 'PA': 'P. A. (P. aeruginosa)'}


def conf(y, yh):
    p = yh >= 0.5
    tp = int((p & (y == 1)).sum()); fn = int((~p & (y == 1)).sum()); fp = int((p & (y == 0)).sum()); tn = int((~p & (y == 0)).sum())
    return metrics(tp, fn, fp, tn)


def metrics(tp, fn, fp, tn):
    sens = tp / (tp + fn) if tp + fn else 0.0; spec = tn / (tn + fp) if tn + fp else 0.0
    den = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return dict(TP=tp, FN=fn, FP=fp, TN=tn, sens=sens, spec=spec, bal_acc=(sens + spec) / 2,
                mcc=float((tp * tn - fp * fn) / den) if den else 0.0)


def q2(y, yh):
    return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)


def rf_oob(Xm, y):
    rf = RandomForestClassifier(300, min_samples_leaf=1, max_features='sqrt', oob_score=True, n_jobs=1, random_state=0).fit(Xm, y)
    return rf.oob_decision_function_[:, 1]


D0 = ng.build('Cleaned_Steps1-4', ycol=YCOL['HI']); X = D0['X']
G, _ = ao.build_groups(X); REPS = np.asarray(ur.centrality_reps(X, G)); assert len(G) == 182
SETS = {'RouteA_all2208': X, 'RouteA_182_centrality': X[:, REPS], 'RouteB_12': core.X_of(core.RED12), 'Morgan2_bits': MB}
out = {}
for k in KEYS:
    y = np.array(P[k]['y'])
    D = ng.build('Cleaned_Steps1-4', ycol=YCOL[k]); assert np.array_equal(D['Y'], y)
    o = {}
    for nm, src in (('PLS-DA, all 2,208 descriptors', 'NoSel2208'), ('PLS-DA, 182 label-free representatives', 'Reps182')):
        o[nm] = conf(y, np.array(P[k]['pred'][src]))
        a = NF['all2208' if src == 'NoSel2208' else 'reps182_centrality'][k]
        assert [o[nm][c] for c in ('TP', 'FN', 'FP', 'TN')] == [a[c] for c in ('TP', 'FN', 'FP', 'TN')], (k, nm)
    for nm, key in (('Random forest, all 2,208 descriptors', 'RouteA_all2208'), ('Random forest, 182 label-free representatives', 'RouteA_182_centrality'),
                    ('Random forest, Route B 12 descriptors', 'RouteB_12'), ('Random forest, Morgan bits', 'Morgan2_bits')):
        po = rf_oob(SETS[key], y)
        assert abs(q2(y, po) - DRF[key][k]['oob_q2']) < 1e-12 and abs(roc_auc_score(y, po) - DRF[key][k]['oob_auc']) < 1e-12, (k, key)
        if key == 'RouteA_all2208':
            assert np.allclose(po, np.array(P[k]['pred']['RF2208']), atol=1e-12), k
        o[nm] = conf(y, po); o[nm]['oob'] = True
    for nm, src in (('PLS-DA, Morgan bits, no selection', 'no_selection'), ('PLS-DA, Morgan bits, filter nested', 'nested_filter')):
        m = MF[k][src]; o[nm] = metrics(m['TP'], m['FN'], m['FP'], m['TN'])
    out[k] = o
    print(k, {nm: (v['TP'], v['FN'], v['FP'], v['TN'], round(v['mcc'], 3)) for nm, v in o.items()}, flush=True)
json.dump(dict(threshold='prediction >= 0.5 counts as active', script='revision2/analysis/classification_controls.py', results=out),
          open(_osp.path.join(_OUT, 'classification_controls.json'), 'w'), indent=1)
print('written', _osp.path.join(_OUT, 'classification_controls.json'))
