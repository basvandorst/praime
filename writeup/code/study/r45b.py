"""Round 45 addendum (derived after the pre-registered formula missed at y = 2, 3, 5): the pre-registered law ignores that
p itself is prime. Given p prime and the first y-rough distance h1 = h1(p mod P_y), the chance that p + h1 is prime is
  (1/lam) prod_{q<=y} (1 - 1/q)^-1 * S_{>y}({0, h1}),  S_{>y}({0,h}) = prod_{q>y} (1 - nu_q/q)/(1 - 1/q)^2
(nu_q = 1 if q | h, else 2), so hit(y) = <S_{>y}({0, h1})> / (lam prod_{q<=y}(1 - 1/q) e^gamma omega(u)).
The average is over the classes of p mod P_y = prod_{q<=y} q coprime to P_y: exact for y <= 13, Monte Carlo (2e5
random classes) for 13 < y <= 2000, and G_y = prod_{q>y}(1 - 2/q)/(1 - 1/q)^2 beyond (where h1 > y is negligible).
usage: r45b.py"""
import os, json, math, random
import numpy as np
import common as C
import theory_r3 as R3
from r45 import buchstab, GAMMA

PR = C.small_primes(20_000_000).astype(np.int64)
LOGG = np.log((1 - 2 / PR[1:].astype(float)) / (1 - 1 / PR[1:].astype(float)) ** 2)      # q >= 3
CUM = np.concatenate([[0.0], np.cumsum(LOGG)])
TOTAL = CUM[-1]


def G(y):
    """prod_{q > y} (1 - 2/q)/(1 - 1/q)^2 (q odd)."""
    k = int(np.searchsorted(PR[1:], y, side='right'))
    return math.exp(TOTAL - CUM[k])


def S_gt(y, h):
    f = G(y)
    for q in PR[1:]:
        q = int(q)
        if q > h:
            break
        if q > y and h % q == 0:
            f *= (1 - 1 / q) / (1 - 2 / q)
    return f


def h1(a, qs):
    h = 1
    while True:
        if all((a + h) % q for q in qs):
            return h
        h += 1


def mean_S(y, rng):
    qs = [int(q) for q in PR if q <= y]
    if y <= 13:
        P = math.prod(qs)
        vals = [S_gt(y, h1(a, qs)) for a in range(1, P) if math.gcd(a, P) == 1]
        return float(np.mean(vals)), 'exact'
    if y <= 2000:
        vals = []
        for _ in range(200000):
            a = [rng.randrange(1, q) for q in qs]               # p mod q, nonzero (p prime > y)
            h = 1
            while True:
                if all((ai + h) % q for ai, q in zip(a, qs)):
                    break
                h += 1
            vals.append(S_gt(y, h))
        return float(np.mean(vals)), 'monte carlo'
    return G(y), 'product'


if __name__ == '__main__':
    r = json.load(open(os.path.join(R3.OUT, 'r45.json')))
    rng = random.Random(45)
    cache = {}
    out = {}
    for name, d in r.items():
        lam = d['lambda']; rows = []
        for row in d['rows']:
            y = row['y']
            if y not in cache:
                cache[y] = mean_S(y, rng)
            ms, how = cache[y]
            pred = min(1.0, row['predicted'] * ms) if row['u'] > 2.0 else 1.0
            rows.append({'y': y, 'u': row['u'], 'measured': row['measured'], 'pre_registered': row['predicted'],
                         'corrected': pred, 'mean_S_gt_y': ms, 'how': how,
                         'diff': row['measured'] - pred, 'sd': row['sd']})
            print(name, 'y %7d meas %.4f prereg %.4f corrected %.4f (%s) diff %+.4f' % (y, row['measured'], row['predicted'], pred, how, row['measured'] - pred), flush=True)
        out[name] = rows
    json.dump(out, open(os.path.join(R3.OUT, 'r45b.json'), 'w'), indent=1)
