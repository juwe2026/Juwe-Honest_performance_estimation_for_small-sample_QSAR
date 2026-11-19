"""MAIN ANALYSIS following the supplementary specification of 24.10.2026, 17:20 (SHA-256 36668775daafd821c1871a15c6a697ea5e3c25afebbb3d6b6a4db818c902df87):
identical with auswertung.py except for datasets_k/, ergebnisse_k/ and bootstrap seed 208 instead of 8.
Originally: evaluation of the ChEMBL benchmark following the frozen pre-specification
(SHA-256 a13bb2c0b513289b160de744bd8d6bc47662d46553ab5a2fe3db6d336a0ae419).

Conventions as in Revision 1 of the manuscript (code archive, Additional file 10):
- permutation p-value  p = (1 + #{null >= observed - 1e-9}) / (1 + N)          (revision/pval_tie.py)
- inflation factor     SD(null, ddof=1) / binomial SD, nominal at 5 % and calibrated at the
                       own null hit rate q = mean(null) / number of tests   (figures/fig7_figS6_hitcount.py)
- paired bootstrap     fixed LOO/OOB predictions, resampled with replacement within each class, 20,000 draws, the same
                       draws for all models of a data set, 95 % percentile interval (revision/paired_bootstrap.py)
                       seed following the supplementary specification: default_rng([93, data set number, 208]).
Number of tests for the hit count: A2 all cleaned descriptors p, A1 and A1u the number of groups G.

Usage: python3 auswertung_k.py  ->  ergebnisse_k/auswertung.json
"""
import json, glob
import numpy as np
from scipy.stats import wilcoxon, spearmanr, page_trend_test, binom

ERG = 'ergebnisse_k'
B = 20000
TIE = 1e-9
BS_CODE = 208   # supplementary specification: bootstrap analysis code 8 + 200
PRIMAER = [f'DS{i:02d}' for i in range(1, 17)]
TEIL = ['DS21', 'DS22', 'DS23', 'DS24']
ROUTEN = ['A1', 'A1u', 'A2']
MODELLE = ['app_A1', 'app_A1u', 'app_A2', 'A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter', 'RF']


def pval(null, obs):
    null = np.asarray(null, float)
    return float((1 + np.sum(null >= obs - TIE)) / (1 + len(null)))


def q2(y, yh):
    return 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)


def q2_matrix(y, yh, I):
    """Q2 for every bootstrap draw (rows of I)."""
    Y = y[I]; P = yh[I]
    return 1 - ((Y - P) ** 2).sum(1) / ((Y - Y.mean(1, keepdims=True)) ** 2).sum(1)


def ziehungen(y, dsnum):
    rng = np.random.default_rng([93, dsnum, BS_CODE])
    i1 = np.where(y == 1)[0]; i0 = np.where(y == 0)[0]
    return np.array([np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]) for _ in range(B)])


def klassifikation(y, p, s=0.5):
    yh = (p >= s).astype(int)
    tp = int(((yh == 1) & (y == 1)).sum()); tn = int(((yh == 0) & (y == 0)).sum())
    fp = int(((yh == 1) & (y == 0)).sum()); fn = int(((yh == 0) & (y == 1)).sum())
    sens = tp / (tp + fn); spez = tn / (tn + fp)
    nenner = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = (tp * tn - fp * fn) / nenner if nenner > 0 else 0.0
    return dict(sens=sens, spez=spez, ba=(sens + spez) / 2, mcc=float(mcc))


def lade_null():
    null = {}
    for f in glob.glob(f'{ERG}/perm/*.json'):
        d = json.load(open(f))
        null.setdefault((d['ds'], d['art']), []).append((d['start'], d['werte']))
    out = {}
    for k, v in null.items():
        v.sort()
        out[k] = [w for _, ws in v for w in ws]
    return out


def main():
    beob = json.load(open(f'{ERG}/beobachtet.json'))
    vorh = json.load(open(f'{ERG}/vorhersagen.json'))
    null = lade_null()
    info = {ds: json.load(open(f'datasets_k/{ds}_info.json'))['info'] for ds in PRIMAER + TEIL}

    tab = {}
    for ds in PRIMAER + TEIL:
        I_ = info[ds]; dsnum = int(ds[2:])
        y = np.array(vorh[ds]['y']); P = {m: np.array(v) for m, v in vorh[ds]['pred'].items()}
        z = dict(zelle=I_['zelle'], n=I_['n'], aktiv=I_['aktiv'], nicht_aktiv=I_['nicht_aktiv'],
                 p=I_['p_sauber'], G=I_['gruppen'], Gp=I_['redundanz_G_p'])
        # models: Q2, AUC, permutation p, classification metrics
        for m in MODELLE:
            q, a = beob[ds][m]
            nl = np.array(null[(ds, m)])
            z[m] = dict(q2=q, auc=a, p_q2=pval(nl[:, 0], q), p_auc=pval(nl[:, 1], a), n_perm=len(nl),
                        **klassifikation(y, P[m]))
        # hit count and inflation
        tests = dict(A1=I_['gruppen'], A1u=I_['gruppen'], A2=I_['p_sauber'])
        nt = null[(ds, 'Treffer')]
        for r in ROUTEN:
            nl = np.array([w[r] for w in nt], float); m_ = tests[r]
            sd = np.std(nl, ddof=1); qn = nl.mean() / m_
            obs = beob[ds]['Treffer'][r]
            z['T_' + r] = dict(beobachtet=obs, null_mittel=float(nl.mean()), null_sd=float(sd), tests=m_,
                               null_rate=float(qn), p=pval(nl, obs),
                               inf_nominal=float(sd / np.sqrt(m_ * 0.05 * 0.95)),
                               inf_kalibriert=float(sd / np.sqrt(m_ * qn * (1 - qn))) if 0 < qn < 1 else float('nan'),
                               perzentil_97_5=float(np.percentile(nl, 97.5)))
        # paired bootstrap
        I = ziehungen(y, dsnum)
        qb = {m: q2_matrix(y, P[m], I) for m in ['A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter', 'RF']}
        bs = {}
        for kontrolle in ['ohneFilter', 'ohneAuswahl', 'RF']:
            for r in ROUTEN:
                d = qb[r] - qb[kontrolle]; lo, hi = np.percentile(d, [2.5, 97.5])
                bs[f'{r}-{kontrolle}'] = dict(diff=float(q2(y, P[r]) - q2(y, P[kontrolle])), lo=float(lo), hi=float(hi),
                                              ganz_ueber_0=bool(lo > 0))
        z['bootstrap'] = bs
        tab[ds] = z
        print(ds, 'done', flush=True)

    # hypotheses on the 16 natural series
    T = {ds: tab[ds] for ds in PRIMAER}
    H = {}
    opt = {r: np.array([T[d]['app_' + r]['q2'] - T[d][r]['q2'] for d in PRIMAER]) for r in ROUTEN}
    H['H1'] = dict(test='Wilcoxon, einseitig, app_A1 - A1 (Q2) > 0',
                   p=float(wilcoxon(opt['A1'], alternative='greater').pvalue),
                   median_diff=float(np.median(opt['A1'])), positiv=int((opt['A1'] > 0).sum()),
                   sekundaer={r: dict(p=float(wilcoxon(opt[r], alternative='greater').pvalue),
                                      median_diff=float(np.median(opt[r])), positiv=int((opt[r] > 0).sum()))
                              for r in ['A1u', 'A2']})
    IF = {r: np.array([T[d]['T_' + r]['inf_kalibriert'] for d in PRIMAER]) for r in ROUTEN}
    pt = page_trend_test(np.column_stack([IF['A1u'], IF['A1'], IF['A2']]))
    H['H2'] = dict(test='Page-Trendtest, Bloecke = 16 Datensaetze, vorhergesagt A1u < A1 < A2 (kalibrierter Faktor)',
                   L=float(pt.statistic), p=float(pt.pvalue), methode=pt.method,
                   sekundaer=dict(A2_gt_A1=float(wilcoxon(IF['A2'] - IF['A1'], alternative='greater').pvalue),
                                  A1_gt_A1u=float(wilcoxon(IF['A1'] - IF['A1u'], alternative='greater').pvalue),
                                  ordnung_vollstaendig=int(((IF['A2'] > IF['A1']) & (IF['A1'] > IF['A1u'])).sum())),
                   median_IF={r: float(np.median(IF[r])) for r in ROUTEN})
    Gp = np.array([T[d]['Gp'] for d in PRIMAER])
    sp = {r: spearmanr(Gp, np.log(IF[r]), alternative='less') for r in ROUTEN}
    H['H3'] = dict(test='Spearman G/p gegen log(kalibrierter Faktor A2), einseitig negativ',
                   rho=float(sp['A2'].statistic), p=float(sp['A2'].pvalue),
                   sekundaer={r: dict(rho=float(sp[r].statistic), p=float(sp[r].pvalue)) for r in ['A1', 'A1u']})
    zaehl = {}
    for kontrolle in ['ohneFilter', 'ohneAuswahl', 'RF']:
        zaehl[kontrolle] = [d for d in PRIMAER if any(T[d]['bootstrap'][f'{r}-{kontrolle}']['ganz_ueber_0'] for r in ROUTEN)]
    k = len(zaehl['ohneFilter'])
    beste = np.array([max(T[d][r]['q2'] for r in ROUTEN) - T[d]['ohneFilter']['q2'] for d in PRIMAER])
    rng = np.random.default_rng([93, 0, BS_CODE])
    med_bs = np.array([np.median(beste[rng.integers(0, 16, 16)]) for _ in range(B)])
    H['H4'] = dict(kriterium='widerlegt, wenn in >= 4 von 16 Datensaetzen ein Intervall (Route - ohneFilter) ganz ueber 0 liegt',
                   anzahl=k, datensaetze=zaehl['ohneFilter'], widerlegt=bool(k >= 4),
                   binomial_p_ab_k=float(binom.sf(k - 1, 16, 0.075)) if k > 0 else 1.0,
                   sekundaer=dict(ohneAuswahl=zaehl['ohneAuswahl'], RF=zaehl['RF']),
                   median_beste_minus_kontrolle=float(np.median(beste)),
                   median_ci95=[float(x) for x in np.percentile(med_bs, [2.5, 97.5])],
                   hinweis_seed=f'Bootstrap ueber Datensaetze: Seed nicht spezifiziert, verwendet default_rng([93, 0, {BS_CODE}])')
    # Holm over H1 to H3 (H4 has its own counting criterion, see report)
    ps = {h: H[h]['p'] for h in ['H1', 'H2', 'H3']}
    for m_ in (3, 4):
        o = sorted(ps, key=ps.get); ok = True; res = {}
        for j, h in enumerate(o):
            ok = ok and ps[h] <= 0.05 / (m_ - j); res[h] = ok
        H[f'holm_m{m_}'] = res
    json.dump(dict(tabelle=tab, hypothesen=H), open(f'{ERG}/auswertung.json', 'w'), indent=1, default=float)
    print(json.dumps(H, indent=1, default=float))


if __name__ == '__main__':
    main()
