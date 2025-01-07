import numpy as np
from scipy import stats

def anderson_ksamp_midrank(samples, Z, Zstar, k, n, N):
    A2akN = 0.
    Z_ssorted_left = Z.searchsorted(Zstar, 'left')
    if N == Zstar.size:
        lj = 1.
    else:
        lj = Z.searchsorted(Zstar, 'right') - Z_ssorted_left
    Bj = Z_ssorted_left + lj / 2.
    for i in np.arange(0, k):
        s = np.sort(samples[i])
        s_ssorted_right = s.searchsorted(Zstar, side='right')
        Mij = s_ssorted_right.astype(float)
        fij = s_ssorted_right - s.searchsorted(Zstar, 'left')
        Mij -= fij / 2.
        inner = lj / float(N) * (N*Mij - Bj*n[i])**2 / (Bj*(N - Bj) - N*lj/4.)
        A2akN += inner.sum() / n[i]
    A2akN *= (N - 1.) / N
    return A2akN


def anderson_ksamp(samples, midrank=True, *, method=None):
    k = len(samples)
    samples = list(map(np.asarray, samples))
    Z = np.sort(np.hstack(samples))
    N = Z.size
    Zstar = np.unique(Z)
    n = np.array([sample.size for sample in samples])

    if midrank:
        A2kN_fun = anderson_ksamp_midrank
    else:
        A2kN_fun = anderson_ksamp_right
    A2kN = A2kN_fun(samples, Z, Zstar, k, n, N)

    return A2kN


df0 = np.array([
    [38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0],
    [39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8],
    [34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0],
    [34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8]
])

res = anderson_ksamp(df0)
print(res)