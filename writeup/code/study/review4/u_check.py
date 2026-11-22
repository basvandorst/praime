"""Own order-3 Bahadur rate for gap g and U_g = int_{1e7}^{P(g)} share_g(x)/ln x dx, for a few g; compare with r24.json."""
import math, sys, json, itertools, time
import numpy as np
from sympy import primerange
sys.path.insert(0, '/home/bas/scratch/study-review4')
from hl_exact import tail_T


def S_rows(rows, g, T):
    rows = np.asarray(rows, dtype=np.int64)
    n, k = rows.shape
    logS = np.zeros(n); dead = np.zeros(n, bool)
    for p in primerange(2, g + 1):
        r = np.sort(rows % p, axis=1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        f = 1 - nu / p
        dead |= f <= 0
        logS += np.log(np.where(f > 0, f, 1.0)) - k * math.log(1 - 1 / p)
    out = np.exp(logS) * T[k]
    out[dead] = 0
    return out


def rate_fn(g):
    T = {k: float(tail_T(k, g)) for k in (2, 3, 4, 5)}
    SB = float(S_rows([[0, g]], g, T)[0])
    I = np.arange(1, g)
    s1 = S_rows(np.stack([np.zeros_like(I), np.full_like(I, g), I], 1), g, T) / SB
    live = I[s1 > 0]; rho = s1[s1 > 0]; n = len(live)
    i, j = np.triu_indices(n, 1)
    s2v = S_rows(np.stack([np.zeros_like(i), np.full_like(i, g), live[i], live[j]], 1), g, T) / SB
    s2 = np.zeros((n, n)); s2[i, j] = s2v; s2[j, i] = s2v
    trip = np.array(list(itertools.combinations(range(n), 3)), dtype=np.int64)
    a, b, c = trip.T
    s3v = np.concatenate([S_rows(np.stack([np.zeros(len(ch), np.int64), np.full(len(ch), g), live[ch[:, 0]], live[ch[:, 1]], live[ch[:, 2]]], 1), g, T)
                          for ch in np.array_split(trip, max(1, len(trip) // 200000))]) / SB

    def L3(lam):
        p = rho / lam
        out = np.sum(np.log1p(-p))
        cc = (s2v / lam ** 2 - p[i] * p[j]) / ((1 - p[i]) * (1 - p[j]))
        out += cc.sum()
        m2 = s2 / lam ** 2
        k3 = s3v / lam ** 3 - m2[a, b] * p[c] - m2[a, c] * p[b] - m2[b, c] * p[a] + 2 * p[a] * p[b] * p[c]
        out -= np.sum(k3 / ((1 - p[a]) * (1 - p[b]) * (1 - p[c])))
        return out
    return SB, n, L3


r24 = json.load(open('/home/bas/praime/writeup/data/study/r24.json'))
Pg = dict(zip(r24['g'], r24['P'])); Ug = dict(zip(r24['g'], r24['U']))
res = {}
for g in [int(x) for x in sys.argv[1:]]:
    t0 = time.time()
    SB, n, L3 = rate_fn(g)
    lg = np.linspace(math.log(1e7), math.log(Pg[g]), 41)
    Lv = np.array([L3(l) for l in lg])
    u = np.linspace(math.log(1e7), math.log(Pg[g]), 4001)
    sh = SB / u * np.exp(np.interp(u, lg, Lv))
    U = float(np.trapezoid(sh * np.exp(u) / u, u))
    res[g] = {'U_own': U, 'U_repo': Ug[g], 'live': n}
    print(g, 'live', n, 'U own %.4f repo %.4f  (%.1fs)' % (U, Ug[g], time.time() - t0), flush=True)
json.dump(res, open('/home/bas/scratch/study-review4/u_check_%s.json' % '_'.join(sys.argv[1:]), 'w'), indent=1)
