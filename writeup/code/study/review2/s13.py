"""Section 12 check, law-free: are the per-residue gap-class tables under-dispersed?
Train window [1e7, 1e8), Q = 210, class c = min(g/2, 20). Split the primes into two halves by the parity of
floor((p - L)/s) (interleaved chunks of length s, so both halves see the same range of ln x), and run a per-row
two-sample chi^2 (homogeneity) test; dof = sum over rows of (cells - 1). Calibration: random halves (must give 1)."""
import math, json, sys
import numpy as np
from sieve import primes

L, W, Q = 10 ** 7, 9 * 10 ** 7, 210
P = primes(L, W)
P = np.append(P, primes(L + W, 2000)[:1])        # the next prime, so that the last gap is defined
g = np.diff(P); p = P[:-1]
cl = np.minimum(g // 2, 20) - 1
units = [a for a in range(Q) if math.gcd(a, Q) == 1]
idx = -np.ones(Q, dtype=int); idx[units] = np.arange(len(units))
st = idx[p % Q]
mask = st >= 0                                     # drop 11, 13 ... none here; p >= 1e7 are all units
st, cl, p = st[mask], cl[mask], p[mask]


def table(sel):
    N = np.zeros((len(units), 20)); np.add.at(N, (st[sel], cl[sel]), 1); return N


def homog(N1, N2):
    chi = 0.0; dof = 0
    for i in range(N1.shape[0]):
        n1, n2 = N1[i].sum(), N2[i].sum(); tot = N1[i] + N2[i]
        E1 = tot * n1 / (n1 + n2); E2 = tot * n2 / (n1 + n2)
        big = (E1 >= 5) & (E2 >= 5)
        o1 = list(N1[i][big]); o2 = list(N2[i][big]); e1 = list(E1[i * 0:][big]) if False else list(E1[big]); e2 = list(E2[big])
        if (~big).any() and E1[~big].sum() >= 5 and E2[~big].sum() >= 5:
            o1.append(N1[i][~big].sum()); o2.append(N2[i][~big].sum()); e1.append(E1[~big].sum()); e2.append(E2[~big].sum())
        o1, o2, e1, e2 = map(np.array, (o1, o2, e1, e2))
        chi += float(np.sum((o1 - e1) ** 2 / e1) + np.sum((o2 - e2) ** 2 / e2)); dof += len(o1) - 1
    return chi / dof, dof


out = {'primes': int(len(p))}
rng = np.random.default_rng(1)
r = [homog(table(s), table(~s))[0] for s in (rng.random(len(p)) < 0.5 for _ in range(10))]
out['random_halves'] = (float(np.mean(r)), float(np.std(r, ddof=1)))
print('random halves', out['random_halves'], flush=True)
for s in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 9 * 10 ** 6, 45 * 10 ** 6):
    sel = ((p - L) // s) % 2 == 0
    v, d = homog(table(sel), table(~sel))
    out['chunks_%d' % s] = (v, d)
    print('interleaved chunks', s, round(v, 3), d, flush=True)
json.dump(out, open('/home/bas/scratch/study-review2/s13.json', 'w'), indent=1)
