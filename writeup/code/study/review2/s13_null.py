"""Calibration of the row-conditioned Pearson chi^2 (pooling E < 5 per row, dof = cells - rows) with multinomial
draws, mod 30 and mod 210 at train, P = empirical row proportions (a stand-in for the law)."""
import math, json
import numpy as np
from sieve import primes

L, W = 10 ** 7, 9 * 10 ** 7
P = primes(L, W + 2000); P = P[: np.searchsorted(P, L + W) + 1]
g = np.diff(P); p = P[:-1]; cl = np.minimum(g // 2, 20) - 1


def chi2(N, Pr):
    E = N.sum(1, keepdims=True) * Pr; chi = 0.0; cells = 0
    for i in range(N.shape[0]):
        big = E[i] >= 5; o = list(N[i][big]); e = list(E[i][big])
        if E[i][~big].sum() >= 5:
            o.append(N[i][~big].sum()); e.append(E[i][~big].sum())
        o = np.array(o); e = np.array(e); e = e * o.sum() / e.sum()
        chi += float(np.sum((o - e) ** 2 / e)); cells += len(o)
    return chi / (cells - N.shape[0]), cells - N.shape[0]


rng = np.random.default_rng(5); out = {}
for Q in (30, 210):
    U = [a for a in range(Q) if math.gcd(a, Q) == 1]; idx = -np.ones(Q, int); idx[U] = np.arange(len(U))
    st = idx[p % Q]; N = np.zeros((len(U), 20)); np.add.at(N, (st, cl), 1)
    Pr = N / N.sum(1, keepdims=True); rows = N.sum(1).astype(int)
    v = [chi2(np.array([rng.multinomial(n, Pr[i]) for i, n in enumerate(rows)], float), Pr)[0] for _ in range(200)]
    out[Q] = {'mean': float(np.mean(v)), 'sd': float(np.std(v, ddof=1)), 'dof': chi2(N, Pr)[1],
              'expected_sd': math.sqrt(2 / chi2(N, Pr)[1])}
    print(Q, out[Q], flush=True)
json.dump(out, open('/home/bas/scratch/study-review2/s13_null.json', 'w'), indent=1)
