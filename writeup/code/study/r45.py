"""Round 45 (preregister/round45.txt): how predictable is the next gap if p mod q is known for every prime q <= y?
Best rule: guess the first y-rough number after p. It is right iff every composite between p and the next prime has
smallest prime factor <= y, so hit(y) = P(M(p) <= y) with M(p) = max spf over the composites in the gap.
Prediction: hit(y) = 1 / (ln x * prod_{q<=y}(1 - 1/q) * e^gamma * omega(u)), u = ln x / ln y (Buchstab).
usage: r45.py"""
import os, json, math
import numpy as np
import common as C
import theory_r3 as R3

GAMMA = 0.5772156649015329


def buchstab(us):
    """omega(u) on a grid: omega = 1/u on [1, 2]; u omega(u) = 1 + int_2^u omega(t - 1) dt beyond."""
    h = 1e-4
    grid = np.arange(1.0, max(us) + 1.0 + h, h)
    om = np.zeros(len(grid))
    i2 = int(round(1.0 / h))
    om[:i2 + 1] = 1 / grid[:i2 + 1]
    F = 1.0
    for i in range(i2 + 1, len(grid)):
        j = i - i2                                           # omega(t - 1) at t = grid[i]
        F += h * 0.5 * (om[j - 1] + om[j])
        om[i] = F / grid[i]
    return np.interp(us, grid, om)


def spf_window(L, W):
    """smallest prime factor of L..L+W-1 for factors up to sqrt(L+W); 0 marks a prime."""
    top = int(math.isqrt(L + W)) + 1
    dt = np.uint16 if top < 65536 else np.uint32
    spf = np.zeros(W, dtype=dt)
    for q in C.small_primes(top):
        q = int(q)
        st = (-L) % q
        s = spf[st::q]
        s[s == 0] = q
        spf[st::q] = s
    return spf


def run(L, W):
    spf = spf_window(L, W)
    pidx = np.nonzero(spf == 0)[0]
    M = np.maximum.reduceat(spf.astype(np.int64), pidx[:-1] + 1)       # each segment ends at the next prime (spf 0)
    gaps = np.diff(pidx)
    lam = math.log(L + W / 2)
    qs = [int(q) for q in C.small_primes(int(math.isqrt(L + W)) + 1)]
    grid = sorted(set([q for q in qs if q <= 50] + [max([q for q in qs if q <= v]) for v in np.geomspace(50, qs[-1], 30)]))
    prod = {}; acc = 1.0; k = 0
    for q in qs:
        acc *= 1 - 1 / q; prod[q] = acc
    us = np.array([lam / math.log(y) for y in grid])
    om = buchstab(us)
    N = len(M)
    rows = []
    for y, u, w in zip(grid, us, om):
        meas = float(np.mean(M <= y))
        pred = min(1.0, 1 / (lam * prod[y] * math.exp(GAMMA) * w))
        mert = 1 / (lam * prod[y])
        sd = math.sqrt(max(meas * (1 - meas), 1e-12) / N)
        rows.append({'y': y, 'u': float(u), 'measured': meas, 'predicted': pred, 'mertens_only': mert, 'sd': sd,
                     'pass': bool(u < 1.5 or abs(meas - pred) <= 0.01 + 3 * sd)})
    # best guess per class mod Q (for comparison with the first-rough-number rule)
    P = pidx[:-1] + L
    best = {}
    for Q in (6, 30, 210):
        r = P % Q; hit = 0
        for a in np.unique(r):
            v, c = np.unique(gaps[r == a], return_counts=True); hit += c.max()
        best[str(Q)] = hit / N
    return {'window': [L, W], 'primes': int(N), 'lambda': lam, 'rows': rows, 'best_guess_per_class': best}


if __name__ == '__main__':
    out = {}
    for name, (L, W) in (('1e9', (10 ** 9, 5 * 10 ** 7)), ('1e12', (10 ** 12, 10 ** 7))):
        out[name] = run(L, W)
        for r in out[name]['rows']:
            print(name, 'y %7d u %6.2f meas %.4f pred %.4f mertens %.4f %s' % (r['y'], r['u'], r['measured'], r['predicted'], r['mertens_only'], '' if r['pass'] else 'MISS'), flush=True)
        print(name, 'best guess per class', out[name]['best_guess_per_class'], flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r45.json'), 'w'), indent=1)
