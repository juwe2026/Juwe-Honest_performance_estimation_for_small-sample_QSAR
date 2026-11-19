"""SENSITIVITY ANALYSIS MEDIAN RULE (addendum 25.10.2026 21:39): as robustheit_k.py, on ergebnisse_m/ with auswertung_m,
purely descriptive: without decision (no "refuted") and without Holm correction. Seeds via auswertung_m (bootstrap 308).
Template: MAIN ANALYSIS (supplementary specification 24.10.2026 17:20): as robustheit.py, on ergebnisse_k/ with auswertung_k.
Robustness check of hypotheses H1 to H4 against numerical instability (decision EW, 24.10.2026).

The frozen analysis (auswertung.py) stays primary. In addition here:
  R1  the same tests without the series with unstable predictions
      (criterion: at least one PLS prediction > 2 or < -1 in any model, result insensitive
      to the choice of the bound between 1.5 and 2, because all values in between lie in series that are flagged anyway)
  R2  stricter and label-free: without the series in which a held-out compound reaches |z| > 100 on any
      descriptor relative to its training compounds (ergebnisse/maxz_labelfrei.json)
  AUC H1 and H4 on all 16 series with the AUC instead of Q2 (the AUC is pre-specified as a test statistic and
      rank-based, a single unstable prediction shifts only the rank of one compound)
H4 bound for fewer than 16 series: smallest number k with P(X >= k) < 0.05 under Binomial(n, 0.075),
the same derivation as in the specification.
Usage: python3 robustheit_m.py  ->  ergebnisse_m/robustheit.json
"""
import json
import numpy as np
from scipy.stats import wilcoxon, spearmanr, page_trend_test, binom, rankdata
import auswertung_m as aw

PRIM = aw.PRIMAER
ROUTEN = aw.ROUTEN


def entgleist(V, ds):
    return any(((np.array(p) > 2) | (np.array(p) < -1)).any() for m, p in V[ds]['pred'].items() if m != 'RF')


def h4_grenze(n):
    return next(k for k in range(1, n + 1) if binom.sf(k - 1, n, 0.075) < 0.05)


def tests(T, dss):
    n = len(dss); R = dict(n=n, datensaetze=dss)
    opt = {r: np.array([T[d]['app_' + r]['q2'] - T[d][r]['q2'] for d in dss]) for r in ROUTEN}
    R['H1'] = {r: dict(p=float(wilcoxon(opt[r], alternative='greater').pvalue), positiv=int((opt[r] > 0).sum()),
                       median=float(np.median(opt[r]))) for r in ROUTEN}
    IF = {r: np.array([T[d]['T_' + r]['inf_kalibriert'] for d in dss]) for r in ROUTEN}
    pt = page_trend_test(np.column_stack([IF['A1u'], IF['A1'], IF['A2']]))
    R['H2'] = dict(p=float(pt.pvalue), ordnung_vollstaendig=int(((IF['A2'] > IF['A1']) & (IF['A1'] > IF['A1u'])).sum()))
    sp = spearmanr([T[d]['Gp'] for d in dss], np.log(IF['A2']), alternative='less')
    R['H3'] = dict(rho=float(sp.statistic), p=float(sp.pvalue))
    k = [d for d in dss if any(T[d]['bootstrap'][f'{r}-ohneFilter']['ganz_ueber_0'] for r in ROUTEN)]
    g = h4_grenze(n)
    beste = np.array([max(T[d][r]['q2'] for r in ROUTEN) - T[d]['ohneFilter']['q2'] for d in dss])
    R['H4'] = dict(anzahl=len(k), datensaetze=k, grenze_hauptanalyse=g,
                   median_beste_minus_kontrolle=float(np.median(beste)))
    return R


def auc_matrix(y, p, I):
    Y = y[I]; P = p[I]
    R_ = rankdata(P, axis=1)
    n1 = Y.sum(1); n0 = Y.shape[1] - n1
    return ((R_ * Y).sum(1) - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    A = json.load(open('ergebnisse_m/auswertung.json')); T = A['tabelle']
    V = json.load(open('ergebnisse_m/vorhersagen.json'))
    Z = json.load(open('ergebnisse_m/maxz_labelfrei.json'))
    aus_R1 = [d for d in PRIM if entgleist(V, d)]
    aus_R2 = [d for d in PRIM if max(Z[d]) > 100]
    out = dict(kriterium_R1='mindestens eine PLS-Vorhersage > 2 oder < -1', ausgeschlossen_R1=aus_R1,
               kriterium_R2='labelfrei: max |z| einer herausgelassenen Verbindung > 100', ausgeschlossen_R2=aus_R2)
    out['alle16'] = tests(T, PRIM)
    out['R1'] = tests(T, [d for d in PRIM if d not in aus_R1])
    out['R2'] = tests(T, [d for d in PRIM if d not in aus_R2])
    # AUC on all 16 series
    da = {r: np.array([T[d]['app_' + r]['auc'] - T[d][r]['auc'] for d in PRIM]) for r in ROUTEN}
    out['AUC_H1'] = {r: dict(p=float(wilcoxon(da[r], alternative='greater').pvalue), positiv=int((da[r] > 0).sum()),
                             median=float(np.median(da[r]))) for r in ROUTEN}
    h4 = {}
    for d in PRIM:
        y = np.array(V[d]['y']); I = aw.ziehungen(y, int(d[2:]))
        P = {m: np.array(V[d]['pred'][m]) for m in ROUTEN + ['ohneFilter']}
        ab = {m: auc_matrix(y, P[m], I) for m in P}
        for r in ROUTEN:
            dd = ab[r] - ab['ohneFilter']; lo, hi = np.percentile(dd, [2.5, 97.5])
            h4[f'{d}_{r}'] = dict(diff=float(T[d][r]['auc'] - T[d]['ohneFilter']['auc']), lo=float(lo), hi=float(hi))
    k = sorted({key[:4] for key, v in h4.items() if v['lo'] > 0})
    beste = np.array([max(T[d][r]['auc'] for r in ROUTEN) - T[d]['ohneFilter']['auc'] for d in PRIM])
    out['AUC_H4'] = dict(anzahl=len(k), datensaetze=k, grenze_hauptanalyse=4,
                         median_beste_minus_kontrolle=float(np.median(beste)), intervalle=h4)
    # check: AUC of the draw matrix for the original sample equals sklearn
    from sklearn.metrics import roc_auc_score
    y = np.array(V['DS01']['y']); p = np.array(V['DS01']['pred']['A1'])
    assert abs(auc_matrix(y, p, np.arange(len(y))[None, :])[0] - roc_auc_score(y, p)) < 1e-12
    json.dump(out, open('ergebnisse_m/robustheit.json', 'w'), indent=1, default=float)
    for key in ['alle16', 'R1', 'R2']:
        R = out[key]
        print(f"\n== {key}: n = {R['n']}")
        print('  H1 Q2:', {r: (round(v['p'], 5), f"{v['positiv']}/{R['n']}", round(v['median'], 3)) for r, v in R['H1'].items()})
        print('  H2:', R['H2'], ' H3:', {k_: round(v, 4) for k_, v in R['H3'].items()})
        print('  H4:', {k_: v for k_, v in R['H4'].items()})
    print('\nexcluded R1:', aus_R1, ' R2:', aus_R2)
    print('AUC H1:', out['AUC_H1'])
    print('AUC H4:', {k_: v for k_, v in out['AUC_H4'].items() if k_ != 'intervalle'})


if __name__ == '__main__':
    main()
