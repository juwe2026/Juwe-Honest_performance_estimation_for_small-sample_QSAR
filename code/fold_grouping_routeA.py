"""Sensitivitaetsanalyse fuer den 2.208-Deskriptor-Pool von Route A.

Tabelle S1 misst die X-Datenabhaengigkeit der unueberwachten Korrelations-
gruppierung fuer die 25 a priori Deskriptoren von Route B. Die Routen A1 und A1u
fixieren dieselbe Art von Schritt auf einem um zwei Groessenordnungen groesseren
Pool: die Partition der 2.208 bereinigten Deskriptoren in 182 Korrelations-
gruppen wird einmal auf allen 31 Verbindungen bestimmt. Hier wird diese Partition
in jedem Trainingsfold aus den 30 Trainingsverbindungen neu gebildet, mit
identischer Regel, und gegen die globale Fassung gestellt.

Die Gruppierung benutzt keine Labels, sie wird deshalb einmal je Fold berechnet
und fuer alle Staemme und alle Permutationen wiederverwendet. Alles
Label-Abhaengige - Repraesentantenwahl nach Effektgroesse (A1), Filter, PLS-DA -
laeuft wie im Manuskript in jedem Fold neu.

Aufruf:  python3 code/fold_grouping_routeA.py [NPERM]
Ausgabe: results/recalc/fold_grouping_routeA.json (QSAR_RECALC_DIR ueberschreibt das)
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_DATA = _os.path.join(_ROOT, "data")
_RES = _os.environ.get("QSAR_RECALC_DIR", _os.path.join(_ROOT, "results", "recalc"))
_os.makedirs(_RES, exist_ok=True)
# Additional file 4 is not bundled here because of its size; place it in data/ under
# this name, or point QSAR_ROUTEA_XLSX at it.
XLSX = "20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx"

import sys, os, json, pickle, time
import numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score

BIG = _os.environ.get("QSAR_ROUTEA_XLSX", _os.path.join(_DATA, XLSX))
CACHE = _os.path.join(_RES, '_foldgroups_A.pkl')
THR = 0.7

STRAINS = ['H.I.', 'S.A.', 'S.Pneu', 'S.Pyo', 'P.A.']
# Die Gruppierung ist labelfrei und fuer alle fuenf Staemme dieselbe, also ist das
# gepoolte Blatt die richtige Quelle sowohl fuer X als auch - ueber seine eigenen
# fuenf Aktivitaetsspalten - fuer Y. Die frueher hier genannten Blaetter CS1-4_*
# gibt es in der aktuellen Datenmatrix nicht mehr; ihre Nachfolger CS1-5_* sind die
# stammspezifischen Saetze NACH dem Effektgroessenfilter und damit nicht der
# ungefilterte 2.208er-Pool, den diese Analyse braucht (Runde-6-Bericht, 8.3).
SHEET = {s: 'Cleaned_Steps1-4' for s in STRAINS}
YCOL = {'H.I.': 'H. I. (H. influenzae)', 'S.A.': 'S. A. (S. aureus)',
        'S.Pneu': 'S. Pneu (S. pneumoniae)', 'S.Pyo': 'S. pyo (S. pyogenes)',
        'P.A.': 'P. A. (P. aeruginosa)'}
FULL = {'H.I.': 'H. influenzae', 'S.A.': 'S. aureus', 'S.Pneu': 'S. pneumoniae',
        'S.Pyo': 'S. pyogenes', 'P.A.': 'P. aeruginosa'}

print("lade Deskriptormatrix ...")
_xl = pd.ExcelFile(BIG)
_clean = pd.read_excel(_xl, 'Cleaned_Steps1-4')
_first = pd.read_excel(_xl, SHEET['H.I.'])
DESC = [c for c in _clean.columns if '::' in c and c in _first.columns]
X = _first[DESC].apply(pd.to_numeric, errors='coerce').values.astype(float)
N, P = X.shape
print("  %d Verbindungen, %d Deskriptoren" % (N, P))
assert P == 2208, "Deskriptorzahl %d, erwartet 2208" % P

Y = {}
for s in STRAINS:
    h = pd.read_excel(_xl, SHEET[s])
    yc = YCOL[s]
    Y[s] = (h[yc].astype(str).str.strip() == 'active').astype(float).values
    assert len(Y[s]) == N
print("  Aktive je Stamm:", {s: int(Y[s].sum()) for s in STRAINS})


# ------------------------------------------------------------------ Gruppierung
def build_groups(Xs, thresh=THR):
    """Identisch zu code/altorder.py build_groups, nur ohne Zwischenspeicher."""
    n, p = Xs.shape
    R = rankdata(Xs, axis=0)

    def adj(M):
        C = np.corrcoef(M, rowvar=False)
        C = np.nan_to_num(C, nan=0.0)
        np.fill_diagonal(C, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            tt = C * np.sqrt((n - 2) / (1 - C ** 2))
        pv = 2 * tdist.sf(np.abs(tt), n - 2)
        pv[np.abs(C) >= 1] = 0.0
        return (np.abs(C) > thresh) & (pv < 0.05)

    ADJ = adj(Xs) | adj(R)
    np.fill_diagonal(ADJ, False)
    var = np.nan_to_num(Xs.var(axis=0))
    remaining = np.ones(p, bool)
    deg = ADJ.sum(1).astype(int)
    groups = []
    while remaining.any():
        d = np.where(remaining, deg, -1)
        mx = d.max()
        if mx <= 0:
            for i in np.where(remaining)[0]:
                groups.append([int(i)])
            break
        cand = np.where(d == mx)[0]
        hub = int(cand[np.argmax(var[cand])])
        members = np.unique(np.concatenate([[hub], np.where(ADJ[hub] & remaining)[0]])).astype(int)
        groups.append([int(m) for m in members])
        for m in members:
            deg -= ADJ[m].astype(int)
        remaining[members] = False
    return groups, ADJ


def centrality_reps(Xs, groups):
    """Repraesentant je Gruppe nach Zentralitaet, wie code/unsup_rep.py."""
    R = rankdata(Xs, axis=0)
    reps = []
    for g in groups:
        g = np.asarray(g)
        if len(g) == 1:
            reps.append(int(g[0])); continue
        Cp = np.corrcoef(Xs[:, g], rowvar=False)
        Cs = np.corrcoef(R[:, g], rowvar=False)
        C = np.maximum(np.abs(np.nan_to_num(Cp)), np.abs(np.nan_to_num(Cs)))
        np.fill_diagonal(C, np.nan)
        reps.append(int(g[int(np.nanargmax(np.nanmean(C, axis=1)))]))
    return np.array(sorted(reps), dtype=int)


def flatten(groups):
    """Gruppen als (order, starts) fuer eine vektorisierte Repraesentantenwahl."""
    order = np.concatenate([np.asarray(g, dtype=int) for g in groups])
    sizes = np.array([len(g) for g in groups], dtype=int)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    return order, starts, sizes


# --------------------------------------------------------- Fold-Gruppierungen
if os.path.exists(CACHE):
    FG = pickle.load(open(CACHE, 'rb'))
    print("Fold-Gruppierungen aus Zwischenspeicher gelesen")
else:
    t0 = time.time()
    GL, _ = build_groups(X)
    print("globale Gruppierung: %d Gruppen (%.1f s)" % (len(GL), time.time() - t0))
    FG = {'global': GL, 'global_cent': centrality_reps(X, GL), 'folds': [], 'folds_cent': []}
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        g, _ = build_groups(X[tr])
        FG['folds'].append(g)
        FG['folds_cent'].append(centrality_reps(X[tr], g))
        print("  Fold %2d: %d Gruppen (%.1f s)" % (i, len(g), time.time() - t0))
    pickle.dump(FG, open(CACHE, 'wb'))

GL = FG['global']
NG = [len(g) for g in FG['folds']]
print("\nglobale Partition: %d Gruppen; foldweise %d bis %d, Median %d"
      % (len(GL), min(NG), max(NG), int(np.median(NG))))

# Wie stark aendert sich die Partition? Anteil der Deskriptorpaare, die in
# beiden Fassungen gemeinsam gruppiert sind beziehungsweise nicht.
def gid_of(groups, p=P):
    gid = np.empty(p, dtype=np.int32)
    for k, g in enumerate(groups):
        gid[np.asarray(g, dtype=int)] = k
    return gid

gid_g = gid_of(GL)
agree = []
for g in FG['folds']:
    gid_f = gid_of(g)
    # Rand-Index ueber eine Stichprobe von Paaren, exakt waere 2.4 Millionen Paare
    rng = np.random.default_rng(7)
    a = rng.integers(0, P, 200000); b = rng.integers(0, P, 200000)
    m = a != b
    agree.append(float(np.mean((gid_g[a[m]] == gid_g[b[m]]) == (gid_f[a[m]] == gid_f[b[m]]))))
print("Rand-Index der Partition, global gegen foldweise: %.4f bis %.4f (Mittel %.4f)"
      % (min(agree), max(agree), float(np.mean(agree))))

FLAT_G = flatten(GL)
FLAT_F = [flatten(g) for g in FG['folds']]

# --------------------------------------------------------------- Pipeline
RANKS = []
for i in range(N):
    tr = np.ones(N, bool); tr[i] = False
    RANKS.append(rankdata(X[tr], axis=0))


def eff_filter(R, ytr, ntr):
    n1 = int(ytr.sum()); n0 = ntr - n1
    U1 = R[ytr.astype(bool)].sum(0) - n1 * (n1 + 1) / 2.0
    rrb = 2 * U1 / (n1 * n0) - 1.0
    mu = n1 * n0 / 2.0
    sigma = np.sqrt(n1 * n0 * (ntr + 1) / 12.0)
    z = (np.abs(U1 - mu) - 0.5) / sigma
    return np.abs(rrb), 2 * norm.sf(z)


def reps_by_eff(flat, eff):
    """Repraesentant je Gruppe = groesste Effektgroesse, vektorisiert."""
    order, starts, sizes = flat
    e = eff[order]
    out = np.empty(len(starts), dtype=int)
    for k, (s, n) in enumerate(zip(starts, sizes)):
        out[k] = order[s + int(np.argmax(e[s:s + n]))]
    return out


def loo(y, mode, per_fold):
    """mode 'A1' (Repraesentant nach Effektgroesse) oder 'A1u' (nach Zentralitaet)."""
    yh = np.zeros(N)
    for i in range(N):
        tr = np.ones(N, bool); tr[i] = False
        ytr = y[tr]; ntr = N - 1
        n1 = int(ytr.sum())
        if n1 < 1 or n1 == ntr:
            yh[i] = ytr.mean(); continue
        R = RANKS[i]
        eff, pv = eff_filter(R, ytr, ntr)
        if mode == 'A1':
            reps = reps_by_eff(FLAT_F[i] if per_fold else FLAT_G, eff)
        else:
            reps = FG['folds_cent'][i] if per_fold else FG['global_cent']
        cols = reps[(eff[reps] > 0.3) & (pv[reps] <= 0.05)]
        if len(cols) < 1:
            yh[i] = ytr.mean(); continue
        Xtr = X[tr][:, cols]; Xte = X[i:i + 1, cols]
        mu = Xtr.mean(0); sd = Xtr.std(0, ddof=1); sd[sd == 0] = 1.0
        m = PLSRegression(n_components=min(1, len(cols), ntr - 1),
                          scale=False).fit((Xtr - mu) / sd, ytr)
        yh[i] = m.predict((Xte - mu) / sd).ravel()[0]
    return yh


def q2(y, yh):
    return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)


def run(s, mode, per_fold, nperm, seed=1234):
    y = Y[s]
    yh = loo(y, mode, per_fold)
    oq, oa = q2(y, yh), roc_auc_score(y, yh)
    if nperm <= 0:
        return dict(q2=float(oq), auc=float(oa))
    qs = np.empty(nperm); as_ = np.empty(nperm)
    for k in range(nperm):
        yp = np.random.default_rng([seed, k]).permutation(y)
        yhp = loo(yp, mode, per_fold)
        qs[k] = q2(yp, yhp); as_[k] = roc_auc_score(yp, yhp)
    return dict(q2=float(oq), auc=float(oa),
                p_q2=float((1 + np.sum(qs >= oq)) / (1 + nperm)),
                p_auc=float((1 + np.sum(as_ >= oa)) / (1 + nperm)))


NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
out = {'n_groups_global': len(GL), 'n_groups_folds': NG,
       'rand_index': {'min': min(agree), 'max': max(agree), 'mean': float(np.mean(agree))},
       'nperm': NPERM, 'results': {}}

print("\n%-8s %-4s %-28s %-28s %s" % ("Stamm", "Rte", "global", "foldweise", "Differenz"))
for mode in ('A1', 'A1u'):
    for s in STRAINS:
        t0 = time.time()
        g = run(s, mode, False, NPERM)
        f = run(s, mode, True, NPERM)
        key = "%s_%s" % (s, mode)
        out['results'][key] = {'global': g, 'fold': f,
                               'dq2': f['q2'] - g['q2'], 'dauc': f['auc'] - g['auc']}
        gp = (" p=%.4f" % g['p_q2']) if NPERM else ""
        fp = (" p=%.4f" % f['p_q2']) if NPERM else ""
        print("%-8s %-4s Q2 %+.4f AUC %.3f%s   Q2 %+.4f AUC %.3f%s   dQ2 %+.4f dAUC %+.4f  (%.0f s)"
              % (s, mode, g['q2'], g['auc'], gp, f['q2'], f['auc'], fp,
                 f['q2'] - g['q2'], f['auc'] - g['auc'], time.time() - t0))

for mode in ('A1', 'A1u'):
    d = [v['dq2'] for k, v in out['results'].items() if k.endswith('_' + mode)]
    out.setdefault('summary', {})[mode] = {'mean_dq2': float(np.mean(d)),
                                           'max_abs_dq2': float(np.max(np.abs(d)))}
    print("%s: mittlere Differenz Q2 %+.4f, groesste absolute %.4f"
          % (mode, np.mean(d), np.max(np.abs(d))))

json.dump(out, open(_os.path.join(_RES, 'fold_grouping_routeA.json'), 'w'), indent=1)
print("geschrieben: %s" % _os.path.join(_RES, "fold_grouping_routeA.json"))
