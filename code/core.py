"""Shared Route B pipeline, identical to code/run_routeB.py and code/lit_plsda.py.

Reads the bundled data/RouteB_descriptor_data.xlsx (English sheet and column names,
matching the published Additional file 2) by default; override with the QSAR_ROUTEB_XLSX
environment variable to point at a different copy.

Corrected 2026-08-27: SRC previously pointed at an absolute, machine-specific path to a file
not included in this archive, and the sheet/column names (German: "Deskriptoren",
"VIP_reduziert", "Reduktion", "S. A. (Staph. Aureus)") no longer matched the published
Additional file 2 (English: "Descriptors", "VIP_reduced", "Reduction", "S. A. (S. aureus)"),
so every script importing this module failed immediately.
"""
import os
import numpy as np, openpyxl
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
SRC = os.environ.get("QSAR_ROUTEB_XLSX", os.path.join(_HERE, "data", "RouteB_descriptor_data.xlsx"))
_wb = openpyxl.load_workbook(SRC, data_only=True)
_rows = list(_wb["Descriptors"].values)
HDR = list(_rows[0]); DATA = [r for r in _rows[1:] if r[0]]
N = len(DATA)
COMPOUNDS = [r[0] for r in DATA]

def _leading_names(rows):
    """First column of a sheet's data rows, stopping at the first blank row (footnote
    text below the table would otherwise be read as a spurious extra descriptor name)."""
    out = []
    for r in rows:
        if not r or r[0] is None:
            break
        out.append(r[0])
    return out

# The 12 representatives that are actually modelled, read from Additional file 2.
RED12 = _leading_names(list(_wb["VIP_reduced"].values)[1:])

# The superseded 13-descriptor set: the same representatives plus C_Count, which was dropped
# once it turned out to be perfectly collinear with EsterLacton_flag (r = rho = 1.000).
#
# It has to be spelled out here rather than read from the sheet. Until now this module read
# the 13-descriptor set from VIP_reduced and filtered C_Count out of it, but that sheet has
# carried no C_Count row since the reduction was documented in it - the row was replaced by a
# footnote. The filter therefore matched nothing, RED13 came out with 12 entries and was
# identical to RED12, and run_all.py's SETS = {"D13": ..., "D12": ...} compared the
# 12-descriptor model with itself. No published number was affected, because every reported
# Route B result uses the 12-descriptor set, but a reader re-running run_all.py would have
# measured a difference of exactly zero between the two sets and could have concluded that
# the collapse of C_Count and EsterLacton_flag makes no difference - which is the opposite of
# what the archived results/results.json shows (round 7 report, 4.1).
#
# C_Count sits in fourth place, the position it occupies in that archived object, so that
# re-running run_all.py reproduces the delivered file rather than a reordered variant.
RED13 = RED12[:3] + ["C_Count"] + RED12[3:]
assert len(RED12) == 12 and len(RED13) == 13, (len(RED12), len(RED13))
assert set(RED13) - set(RED12) == {"C_Count"}
def _reduction_data_rows(rows):
    """Reduction sheet's data rows only, stopping at the first fully blank row (an
    explanatory Note. paragraph follows the table and must not be read as data)."""
    out = []
    for r in rows:
        if not r or all(v is None for v in r):
            break
        out.append(r)
    return out

_RED_ROWS = _reduction_data_rows(list(_wb["Reduction"].values)[1:])
LIT25 = [r[0] for r in _RED_ROWS if r[0]] + \
        [p.strip().split(" (")[0] for r in _RED_ROWS if r[0] and r[2]
         for p in str(r[2]).split(",")]

STRAINS = ['H.I.', 'S.A.', 'S.Pneu', 'S.Pyo', 'P.A.']
YCOL = {'H.I.':'H. I. (H. influenzae)', 'S.A.':'S. A. (S. aureus)',
        'S.Pneu':'S. Pneu (S. pneumoniae)', 'S.Pyo':'S. pyo (S. pyogenes)',
        'P.A.':'P. A. (P. aeruginosa)'}
FULLNAME = {'H.I.':'H. influenzae', 'S.A.':'S. aureus', 'S.Pneu':'S. pneumoniae',
            'S.Pyo':'S. pyogenes', 'P.A.':'P. aeruginosa'}
# nLV per strain for the published VIP table: uniformly 2 components for every strain
# (matching Fig. S2b/Fig. S3a and Additional file 2's VIP_reduced sheet). Previously this
# used a mixed 4/4/3/1 per strain that matched neither the figures nor the validated (nLV=1)
# models; corrected together with Additional file 2's VIP_reduced sheet.
VIP_NLV = {'H.I.': 2, 'P.A.': 2, 'S.A.': 2, 'S.Pneu': 2, 'S.Pyo': 2}

def col(name):
    j = HDR.index(name); return np.array([float(r[j]) for r in DATA])
def X_of(cols):
    return np.column_stack([col(c) for c in cols])
def y_of(s):
    j = HDR.index(YCOL[s])
    return np.array([1.0 if str(r[j]).strip().lower() == 'active' else 0.0 for r in DATA])

def loo(X, y, nlv):
    yh = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        Xt = X[tr]; mu = Xt.mean(0); sd = Xt.std(0, ddof=1); sd[sd == 0] = 1
        m = PLSRegression(n_components=min(nlv, X.shape[1], N-2), scale=False).fit((Xt-mu)/sd, y[tr])
        yh[i] = m.predict((X[i:i+1]-mu)/sd).ravel()[0]
    return yh

def q2(y, yh):
    return 1 - np.sum((y-yh)**2)/np.sum((y-y.mean())**2)

def autoscale_all(X):
    mu = X.mean(0); sd = X.std(0, ddof=1); sd[sd == 0] = 1; return (X-mu)/sd

def vip(m, X):
    t, w, q = m.x_scores_, m.x_weights_, m.y_loadings_
    p, h = X.shape[1], t.shape[1]
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    return np.sqrt(p * ((w/np.linalg.norm(w, axis=0))**2 @ s).ravel() / s.sum())

def vip_table(cols, nlv_map=None):
    X = autoscale_all(X_of(cols)); out = {}
    for s in STRAINS:
        nlv = (nlv_map or VIP_NLV)[s]
        m = PLSRegression(n_components=nlv, scale=False).fit(X, y_of(s))
        out[s] = dict(zip(cols, [float(v) for v in vip(m, X)]))
    return out

def perm_test(cols, s, nlv, nperm=2000, seed=702):
    """Same seed convention as code/run_routeB.py."""
    X = X_of(cols); y = y_of(s)
    yh = loo(X, y, nlv); oq = q2(y, yh); oa = roc_auc_score(y, yh)
    qs = np.empty(nperm); as_ = np.empty(nperm)
    for k in range(nperm):
        yp = np.random.default_rng([seed, k]).permutation(y)
        yhp = loo(X, yp, nlv); qs[k] = q2(yp, yhp); as_[k] = roc_auc_score(yp, yhp)
    return dict(q2=float(oq), auc=float(oa),
                p_q2=float((1+np.sum(qs >= oq))/(1+nperm)),
                p_auc=float((1+np.sum(as_ >= oa))/(1+nperm)),
                n_ge_q2=int(np.sum(qs >= oq)), n_ge_auc=int(np.sum(as_ >= oa)))

def confusion(cols, s, nlv=1, thr=0.5):
    y = y_of(s); yh = loo(X_of(cols), y, nlv); p = (yh >= thr).astype(float)
    tp = int(((p==1)&(y==1)).sum()); fn = int(((p==0)&(y==1)).sum())
    fp = int(((p==1)&(y==0)).sum()); tn = int(((p==0)&(y==0)).sum())
    den = np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return dict(TP=tp, FN=fn, FP=fp, TN=tn,
                sens=tp/(tp+fn) if tp+fn else 0.0,
                spec=tn/(tn+fp) if tn+fp else 0.0,
                bacc=0.5*((tp/(tp+fn) if tp+fn else 0)+(tn/(tn+fp) if tn+fp else 0)),
                mcc=float((tp*tn-fp*fn)/den) if den else 0.0)
