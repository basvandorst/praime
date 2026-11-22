"""Round 91 (preregister/round91.txt): Montgomery's pair correlation tested with primes only, through Goldston-Montgomery:
R(X, delta) = int_1^X (psi((1+delta)x) - psi(x) - delta x)^2 dx / ((1/2) delta X^2) against ln(1/delta) + B,
B = 1 - gamma - ln 2pi. psi exact (prime powers included) from a segmented sieve; the integral by the midpoint rule.
usage: R91_X=1e9 r91.py"""
import os, json, math, time
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
X = int(float(os.environ.get('R91_X', 1e9)))
DELTAS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6] + ([1e-7] if X >= 10 ** 10 else [])
NPT = 2_000_000
STEP = X / NPT
B = 1 - 0.5772156649015329 - math.log(2 * math.pi)


def psi_at(q):
    """psi at the sorted float array q (values <= qmax), exact up to float rounding: sum of ln p over p^k <= x."""
    qmax = int(q[-1]) + 1
    out = np.zeros(len(q)); base = 0.0; SEG = 10 ** 7
    # prime powers p^k, k >= 2: few, added separately
    small = C.small_primes(int(math.isqrt(qmax)) + 2).astype(np.int64)
    pw = []
    for p in small.tolist():
        v = p * p
        while v <= qmax:
            pw.append((v, math.log(p))); v *= p
    pw.sort(); pwv = np.array([a for a, _ in pw], dtype=np.float64); pwl = np.cumsum([b for _, b in pw]) if pw else np.zeros(0)
    for s0 in range(0, qmax + 1, SEG):
        n = min(SEG, qmax + 1 - s0)
        if s0 == 0:
            m = np.zeros(n, bool); P = C.small_primes(n).astype(np.int64); P = P[P < n]
        else:
            P = np.nonzero(C.primes_mask(s0, n))[0].astype(np.int64) + s0
        lp = np.log(P.astype(np.float64)); cs = np.cumsum(lp)
        lo = np.searchsorted(q, s0, side='left'); hi = np.searchsorted(q, s0 + n, side='left')
        if hi > lo:
            k = np.searchsorted(P, np.floor(q[lo:hi]).astype(np.int64), side='right')
            out[lo:hi] = base + np.concatenate([[0.0], cs])[k]
        base += math.fsum(lp.tolist()) if len(lp) < 2_000_000 else float(np.sum(lp))
    # add prime powers
    if len(pwv):
        k = np.searchsorted(pwv, np.floor(q), side='right')
        out += np.where(k > 0, pwl[np.maximum(k - 1, 0)], 0.0)
    return out


if __name__ == '__main__':
    t0 = time.time()
    STEP_ = (X - 1) / NPT
    # stratified random points (one per cell, seed 91): a grid commensurate with the integers biased small delta
    # (twenty-eighth review: every window started in the same residue class)
    x = (np.arange(NPT) + np.random.default_rng(91).random(NPT)) * STEP_ + 1.0
    xs = x[x <= X]
    pts = [xs] + [xs * (1 + d) for d in DELTAS]
    allq = np.concatenate(pts); order = np.argsort(allq, kind='stable')
    vals = np.empty_like(allq); vals[order] = psi_at(allq[order])
    n = len(xs); ps0 = vals[:n]
    # sanity: psi(X) against X
    res = {'X': X, 'step': STEP, 'B': B, 'points': 'stratified random, seed 91', 'rows': []}
    for i, d in enumerate(DELTAS):
        ps1 = vals[(i + 1) * n:(i + 2) * n]
        diff = ps1 - ps0 - d * xs
        I = float(np.sum(diff * diff) * STEP_)
        R = I / (0.5 * d * X * X)
        res['rows'].append({'delta': d, 'R': R, 'gue_second_order': math.log(1 / d) + B, 'gue_leading': math.log(1 / d),
                            'cramer': math.log(X) - 1.5, 'diff_second_order': R - (math.log(1 / d) + B)})
        print('delta %.0e  R %.4f  ln(1/d)+B %.4f  diff %+.4f  ln(1/d) %.3f  Cramer %.2f' % (d, R, math.log(1 / d) + B, R - math.log(1 / d) - B, math.log(1 / d), math.log(X) - 1.5), flush=True)
    res['psi_at_last_midpoint'] = float(ps0[-1]); res['last_midpoint'] = float(xs[-1]); res['seconds'] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, 'r91_X%.0e.json' % X), 'w'), indent=1)
    print('psi(%.0f) = %.1f, %d s' % (xs[-1], ps0[-1], time.time() - t0))
