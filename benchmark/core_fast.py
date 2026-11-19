"""Vectorised version of the univariate filter.
Computes the rank-biserial effect size and the Mann-Whitney p for ALL descriptors at once.
Must be numerically identical to scipy.stats.mannwhitneyu (asymptotic, with tie and
continuity correction). The check against it is in pruefe_core_fast().
"""
import numpy as np
from scipy.stats import rankdata, norm


def effekt_und_p(X, y):
    """X: (n, p) matrix, y: (n,) with 0/1. Returns (r_rb, p) per column."""
    n = X.shape[0]
    n1 = int(y.sum()); n2 = n - n1
    R = rankdata(X, axis=0)                      # mid-rank, as scipy
    R1 = R[y == 1].sum(axis=0)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    r_rb = 2.0 * U1 / (n1 * n2) - 1.0
    # tie correction per column
    tie_term = np.zeros(X.shape[1])
    for j in range(X.shape[1]):
        _, cnt = np.unique(R[:, j], return_counts=True)
        cnt = cnt[cnt > 1]
        if cnt.size:
            tie_term[j] = (cnt ** 3 - cnt).sum()
    mu = n1 * n2 / 2.0
    sigma = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1.0))))
    U = np.maximum(U1, n1 * n2 - U1)             # two-sided: larger statistic
    with np.errstate(invalid='ignore', divide='ignore'):
        z = (U - mu - 0.5) / sigma
        p = 2.0 * norm.sf(z)
    p = np.clip(np.nan_to_num(p, nan=1.0), 0.0, 1.0)
    return r_rb, p


def pruefe_core_fast(X, y, stichprobe=60, seed=7):
    """Compares the vectorised version column by column with scipy."""
    from scipy.stats import mannwhitneyu
    rng = np.random.default_rng(seed)
    cols = rng.choice(X.shape[1], size=min(stichprobe, X.shape[1]), replace=False)
    r_fast, p_fast = effekt_und_p(X, y)
    n1 = int(y.sum()); n2 = len(y) - n1
    dr = dp = 0.0
    for j in cols:
        u, p = mannwhitneyu(X[y == 1, j], X[y == 0, j], alternative='two-sided', method='asymptotic')
        r = 2 * u / (n1 * n2) - 1
        dr = max(dr, abs(r - r_fast[j]))
        dp = max(dp, abs(p - p_fast[j]))
    return dict(spalten=len(cols), max_abw_effektstaerke=float(dr), max_abw_p=float(dp))
