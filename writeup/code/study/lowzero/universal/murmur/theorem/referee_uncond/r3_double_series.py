"""Referee: E(r) = sum_{m,n} e(m,n)/(mn) by direct partial sums (m <= 30, n <= X), r = 0..3, N = 20011,
e(m,n) from local means computed by enumeration (l | 2m, l | r) and closed forms otherwise; against nu'(r), zeta(2)(1-N^-2).
Tail of the n-sum is ~ c X^(-1/2) (square parts): extrapolate from X = 10^5/4, 10^5 (Richardson in X^(-1/2))."""
import math
import numpy as np
from sympy import factorint, primerange

N = 20011
X = 10 ** 5
spf = np.zeros(X + 1, dtype=np.int64)
for l in primerange(2, X + 1):
    blk = spf[l::l]; blk[blk == 0] = l


def kron2(D):
    return np.where(D % 2 == 0, 0, np.where((D % 8 == 1) | (D % 8 == 7), 1, -1))


def leg(D, l):
    return np.array([0 if x % l == 0 else (1 if pow(int(x), (l - 1) // 2, l) == 1 else -1) for x in D])


cache = {}


def e_loc(l, i, jc, r):
    key = (l, i, jc, r)
    if key in cache:
        return cache[key]
    if l == N:
        v = 1.0 if (i == 0 and jc == 0) else 0.0
    else:
        mod = l ** (2 * i + 1) if l != 2 else 2 ** (2 * i + 4)
        ps = np.array([p for p in range(mod) if p % l], dtype=np.int64)
        ps = ps + mod * ((r * r * N) // (4 * mod) + 2)
        D = N * (4 * ps - r * r * N)
        ok = D % (l ** (2 * i)) == 0
        Dl = -(D // (l ** (2 * i)))
        if l == 2:
            ok &= (Dl % 4 == 0) | (Dl % 4 == 1)
            kap = kron2(Dl)
        else:
            kap = leg(Dl, l)
        w = ok * (1 if jc == 0 else (kap if jc == 'o' else kap * kap))
        v = float(w.mean())
    cache[key] = v
    return v


def E_closed(l, r):
    if l == N:
        return 1.0
    if l == 2:
        return 2 / 3 if r % 2 else 4 / 3
    if r % l == 0:
        return l * l / (l * l - 1)
    return 1 - 1 / ((l - 1) ** 2 * (l + 1))


for r in range(4):
    tot = {X // 4: 0.0, X: 0.0}
    for m in range(1, 31):
        fm = factorint(m)
        if r == 0 and m > 2:
            continue
        em = 1.0
        for l, i in fm.items():
            em *= e_loc(l, i, 0, r)
        if r % 2 == 0 and 2 not in fm:
            em *= 1.0
        if em == 0:
            continue
        # multiplicative ratio g(n) = e(m,n)/e(m,1): local ratio at l^j || n
        g = np.zeros(X + 1); g[1] = 1.0
        for n in range(2, X + 1):
            l = int(spf[n]); j = 0; t = n
            while t % l == 0:
                t //= l; j += 1
            i = fm.get(l, 0)
            jc = 'o' if j % 2 else 'e'
            if l == N:
                loc = 0.0
            elif i == 0 and l != 2:
                if r % l == 0:  # l | r (r = 0 included): kappa = +-1 equally
                    loc = 0.0 if jc == 'o' else 1.0
                else:
                    loc = -1 / (l - 1) if jc == 'o' else 1 - 1 / (l - 1)
            else:
                base = e_loc(l, i, 0, r)
                loc = e_loc(l, i, jc, r) / base if base else 0.0
            g[n] = g[t] * loc
        nn = np.arange(1, X + 1)
        cs = np.cumsum(g[1:] / nn)
        for Xi in tot:
            tot[Xi] += em / m * cs[Xi - 1]
    E = 1.0; nu = 1.0
    for l in primerange(2, 2 * 10 ** 6):
        E *= E_closed(l, r)
        nu *= (l * l / (l * l - 1)) if (r % l == 0) else (1 - 1 / ((l - 1) ** 2 * (l + 1)))
    a, b = tot[X // 4], tot[X]
    extr = 2 * b - a  # tail ~ c X^-1/2: S(X) = S - c X^-1/2, S(X/4) = S - 2c X^-1/2
    ref = math.pi ** 2 / 6 * (1 - N ** -2) if r == 0 else nu
    print(f"r={r}: partial S(X/4)={a:.5f} S(X)={b:.5f} extrapolated {extr:.5f} | E(r) {E:.5f} | "
          f"{'zeta(2)(1-N^-2)' if r == 0 else 'nu(r)'} {ref:.5f}")
