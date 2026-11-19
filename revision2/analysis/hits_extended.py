# -*- coding: utf-8 -*-
"""Monte Carlo precision check of the univariate hit-count p-values (response to Reviewer 2, point 2).

The permutation streams of the archived hit-count analyses are extended from 1,000 to B = 100,000
label permutations. B was fixed at 100,000 for every hit-count test of all three data sets before
this script was run on the peptide and Gaucher data. The schemes are those of the archived runs,
so the first 1,000 permutations are identical to the archived ones (checked below):
  monoterpenoids  code/make_hits_three_variants.py, a new numpy default_rng(9) for each strain,
                  successive calls of permutation(y)
  peptides        revision/lone/run2_tie.py, default_rng([9207, j]) (main definition) and
                  default_rng([9307, j]) (sensitivity definition) for permutation j
  Gaucher         revision/jobG_gaucher_hits.py, default_rng([8801, j])
  tie-break rules revision/jobB_hits.py (Table S6, largest variance and first index), default_rng([8202, j])
For every test the script reports the observed count, the permutation p-value (k + 1)/(B + 1) for
B = 1,000, 10,000 and 100,000, the exceedance count k and the 95 % Clopper-Pearson interval of
k/B. The null mean, the null standard deviation and the calibrated inflation factor are given for
B = 1,000 (the values used in the article) and for B = 100,000 (for comparison only).
Usage: python hits_extended.py {mono|pep|gau|tie|all} [B]   (writes hits_extended_<part>.json to _OUT)
"""
import os as _osp, sys
_ROOT = _osp.path.dirname(_osp.path.abspath(__file__))
while not _osp.path.exists(_osp.path.join(_ROOT, 'CHANGELOG_revision2.md')) and _osp.path.dirname(_ROOT) != _ROOT:
    _ROOT = _osp.path.dirname(_ROOT)
_OUT = _osp.environ.get('QSAR_REV2_OUT', _osp.path.join(_ROOT, 'revision2', 'results'))
sys.dont_write_bytecode = True
sys.path.insert(0, _osp.path.join(_ROOT, 'code'))
import json, time
import numpy as np, pandas as pd
from scipy.stats import rankdata, beta

PART = sys.argv[1] if len(sys.argv) > 1 else 'all'
B = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
CHECK_B = (1000, 10000, 100000)


def cp(k, n):
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return lo, hi


def summarise(obs, null, ntests):
    d = dict(n_tests=int(ntests), obs=int(obs))
    for b in CHECK_B:
        if b > len(null):
            continue
        k = int((null[:b] >= obs).sum()); lo, hi = cp(k, b)
        d['B%d' % b] = dict(k=k, p=(k + 1) / (b + 1), ci95_k_over_B=[lo, hi])
    for b in (1000, len(null)):
        nl = null[:b]; m = float(nl.mean()); s = float(nl.std(ddof=1)); q = m / ntests
        d['null_B%d' % b] = dict(mean=m, sd=s, rate=q, infl_own=s / np.sqrt(ntests * q * (1 - q)) if 0 < q < 1 else None)
    return d


def write(name, res):
    os = _osp
    path = os.path.join(_OUT, 'hits_extended_%s.json' % name)
    json.dump(dict(B=B, script='revision2/analysis/hits_extended.py', results=res), open(path, 'w'), indent=1)
    print('written', path, flush=True)


def mono():
    import nested_generic as ng, altorder as ao, unsup_rep as ur
    import importlib.util
    spec = importlib.util.spec_from_file_location('mh', _osp.path.join(_ROOT, 'code', 'make_hits_three_variants.py'))
    mh = importlib.util.module_from_spec(spec); spec.loader.exec_module(mh)
    ARCH = np.load(_osp.path.join(_ROOT, 'results', 'raw_permutations', 'hits_three_variants.npy'), allow_pickle=True).item()
    D0 = ng.build('Cleaned_Steps1-4', ycol=mh.YCOL['H.I.']); X = D0['X']
    groups, _ = ao.build_groups(X); assert len(groups) == 182
    reps = np.asarray(ur.centrality_reps(X, groups)); R = rankdata(X, axis=0)
    res = {}
    for lab, col in mh.YCOL.items():
        t = time.time()
        D = ng.build('Cleaned_Steps1-4', ycol=col); assert np.array_equal(D['X'], X)
        y = D['Y']; N = D['N']
        obs = mh.hit_counts(R, y, N, groups, reps)
        rng = np.random.default_rng(9)
        null = np.array([mh.hit_counts(R, rng.permutation(y), N, groups, reps) for _ in range(B)], dtype=np.int64)
        for i, k in enumerate(['null_A2', 'null_A1s', 'null_A1u']):
            assert np.array_equal(null[:1000, i], ARCH[lab][k]), (lab, k)
        assert (obs[0], obs[1], obs[2]) == (ARCH[lab]['obs_A2'], ARCH[lab]['obs_A1s'], ARCH[lab]['obs_A1u'])
        res[lab] = {nm: summarise(obs[i], null[:, i], X.shape[1] if nm == 'A2' else len(groups))
                    for i, nm in enumerate(['A2', 'A1', 'A1u'])}
        res[lab]['first_1000_identical_to_archive'] = True
        print(lab, {nm: (res[lab][nm]['obs'], round(res[lab][nm]['B1000']['p'], 4), round(res[lab][nm]['B%d' % B]['p'], 5)) for nm in ['A2', 'A1', 'A1u']},
              '%.0f s' % (time.time() - t), flush=True)
        write('monoterpenoids', res)


def pep():
    import nested_generic as ng, altorder as ao, unsup_rep as ur
    REV = _osp.path.join(_ROOT, 'revision')
    rows = json.load(open(_osp.path.join(REV, 'lone', 'rows.json')))
    Xd = pd.read_csv(_osp.path.join(REV, 'lone', 'desc_pool.csv')); assert list(Xd['No'].astype(str)) == [r['No'] for r in rows]
    DESC = [c for c in Xd.columns if c != 'No']; X = Xd[DESC].values.astype(float); N = X.shape[0]
    G, _ = ao.build_groups(X); REPS = np.asarray(ur.centrality_reps(X, G)); R = rankdata(X, axis=0)
    def hits(yy):
        e, p = ng.eff_filter(R, yy, N); ps = (e > 0.3) & (p <= 0.05)
        a1 = np.array([g[int(np.argmax(e[g]))] for g in G])
        return int(ps.sum()), int(ps[a1].sum()), int(ps[REPS].sum())
    res = {}
    for rule, seed in (('A', 9207), ('C', 9307)):
        THR = {'A': {'EC': 8, 'SA': 8, 'PA': 16, 'AB': 8}, 'C': {'EC': 16, 'SA': 16, 'PA': 16, 'AB': 16}}[rule]
        ARCH = json.load(open(_osp.path.join(_ROOT, 'results', 'revision', 'lone', 'res_%s.json' % rule)))
        res[rule] = {}
        for k in ['EC', 'SA', 'PA', 'AB']:
            t = time.time()
            y = np.array([1.0 if (float(r['MIC_' + k]) <= THR[k] and r['censored_' + k] != 'yes') else 0.0 for r in rows])
            obs = hits(y)
            null = np.array([hits(np.random.default_rng([seed, j]).permutation(y)) for j in range(B)], dtype=np.int64)
            d = {}
            for i, nm in enumerate(['A2', 'A1', 'A1u']):
                d[nm] = summarise(obs[i], null[:, i], X.shape[1] if nm == 'A2' else len(G))
                a = ARCH[k]['hits'][nm]
                assert a['obs'] == obs[i] and abs(a['p'] - d[nm]['B1000']['p']) < 1e-12, (rule, k, nm)
                assert abs(a['null_mean'] - d[nm]['null_B1000']['mean']) < 1e-9 and abs(a['null_sd'] - d[nm]['null_B1000']['sd']) < 1e-9
            d['first_1000_match_archive'] = 'observed count, p, null mean and null sd equal to res_%s.json' % rule
            res[rule][k] = d
            print(rule, k, {nm: (d[nm]['obs'], round(d[nm]['B1000']['p'], 4), round(d[nm]['B%d' % B]['p'], 5)) for nm in ['A2', 'A1', 'A1u']},
                  '%.0f s' % (time.time() - t), flush=True)
            write('peptides', res)


def gau():
    import gaucher_pipe as gp
    X, y, mz, _ = gp.load(_osp.path.join(_ROOT, 'data', 'Gaucherdata_paper39_median_normalised.csv')); X, mz, _ = gp.clean(X, mz)
    n, p = X.shape; G, _ = gp.build_groups(X); reps = np.array(gp.centrality_reps(X, G)); R = rankdata(X, axis=0)
    def counts(yy):
        e, pv = gp.eff_filter(R, yy, n); ps = (e > 0.3) & (pv <= 0.05)
        a1 = np.array([g[int(np.argmax(e[g]))] for g in G])
        return int(ps.sum()), int(ps[a1].sum()), int(ps[reps].sum())
    ARCH = json.load(open(_osp.path.join(_ROOT, 'results', 'revision', 'G_gaucher_hits.json')))
    t = time.time(); obs = counts(y)
    null = np.array([counts(np.random.default_rng([8801, j]).permutation(y)) for j in range(B)], dtype=np.int64)
    res = {}
    for i, (nm, nt) in enumerate([('A2', p), ('A1', len(G)), ('A1u', len(G))]):
        res[nm] = summarise(obs[i], null[:, i], nt)
        a = ARCH[nm]; assert a['obs'] == obs[i] and abs(a['p'] - res[nm]['B1000']['p']) < 1e-12, nm
    res['first_1000_match_archive'] = 'observed count and p equal to G_gaucher_hits.json'
    print('Gaucher', {nm: (res[nm]['obs'], round(res[nm]['B1000']['p'], 4), round(res[nm]['B%d' % B]['p'], 5)) for nm in ['A2', 'A1', 'A1u']},
          '%.0f s' % (time.time() - t), flush=True)
    write('gaucher', res)


def tie():
    sys.path.insert(0, _osp.path.join(_ROOT, 'revision'))
    from common import X, Y, KEYS, N, ADJ
    import nested_generic as ng
    p = X.shape[1]; R = rankdata(X, axis=0)
    Cp = np.abs(np.nan_to_num(np.corrcoef(X, rowvar=False))); Cs = np.abs(np.nan_to_num(np.corrcoef(R, rowvar=False)))
    C = np.maximum(Cp, Cs); np.fill_diagonal(C, np.nan)
    def groups(tb):
        rem = np.ones(p, bool); deg = ADJ.sum(1).astype(int); Gs = []; var = X.var(0)
        while rem.any():
            d = np.where(rem, deg, -1); mx = d.max()
            if mx <= 0:
                Gs += [[int(i)] for i in np.where(rem)[0]]; break
            cand = np.where(d == mx)[0]
            hub = cand[np.argmax(var[cand])] if tb == 'var' else cand[0]
            mem = np.unique(np.concatenate([[hub], np.where(ADJ[hub] & rem)[0]])); Gs.append(mem.tolist())
            for m_ in mem: deg -= ADJ[m_].astype(int)
            rem[mem] = False
        return Gs
    def creps(Gs):
        return np.array([g[0] if len(g) == 1 else g[int(np.nanargmax(np.nanmean(C[np.ix_(g, g)], 1)))] for g in Gs])
    GV = {'variance (as published)': groups('var'), 'first index': groups('first')}
    RP = {g: creps(v) for g, v in GV.items()}
    assert len(GV['variance (as published)']) == 182
    ARCH = json.load(open(_osp.path.join(_ROOT, 'results', 'revision', 'B_hits.json')))['tiebreak_fixed']
    res = {g: {'n_groups': len(GV[g])} for g in GV}
    for k in KEYS:
        t = time.time(); y = Y[k]
        e, pv = ng.eff_filter(R, y, N); po = (e > 0.3) & (pv <= 0.05)
        null = {g: np.empty(B, dtype=np.int64) for g in GV}
        for j in range(B):
            e_, pv_ = ng.eff_filter(R, np.random.default_rng([8202, j]).permutation(y), N); ps = (e_ > 0.3) & (pv_ <= 0.05)
            for g in GV: null[g][j] = int(ps[RP[g]].sum())
        for g in GV:
            o = int(po[RP[g]].sum()); d = summarise(o, null[g], len(GV[g])); a = ARCH[g][k]
            assert a['obs'] == o and abs(a['p'] - d['B1000']['p']) < 1e-12 and abs(a['null_mean'] - d['null_B1000']['mean']) < 1e-9, (g, k)
            res[g][k] = d
        print(k, {g: (res[g][k]['obs'], round(res[g][k]['B1000']['p'], 4), round(res[g][k]['B%d' % B]['p'], 5)) for g in GV}, '%.0f s' % (time.time() - t), flush=True)
        write('tiebreak', res)


if __name__ == '__main__':
    if PART in ('mono', 'all'): mono()
    if PART in ('pep', 'all'): pep()
    if PART in ('gau', 'all'): gau()
    if PART in ('tie', 'all'): tie()
