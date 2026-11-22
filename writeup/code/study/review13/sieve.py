"""Own segmented odd-only sieve; counts gaps >= 200 between consecutive primes in [L0, L1).
usage: sieve.py L0 L1 out.json  (segments of 1e8 numbers)"""
import sys, json, math, time
import numpy as np


def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def odd_primes_in(lo, hi, bp):
    """primes in [lo, hi), lo and hi even, lo > max(bp)^2 not required but lo > 2"""
    n = (hi - lo) // 2            # index i <-> lo + 1 + 2i
    m = np.ones(n, bool)
    for p in bp[1:]:
        p = int(p)
        if p * p >= hi:
            break
        # first odd multiple of p >= max(p*p, lo+1)
        st = max(p * p, ((lo + 1 + p - 1) // p) * p)
        if st % 2 == 0:
            st += p
        m[(st - lo - 1) // 2::p] = False
    return lo + 1 + 2 * np.nonzero(m)[0].astype(np.int64)


if __name__ == '__main__':
    L0, L1, out = int(float(sys.argv[1])), int(float(sys.argv[2])), sys.argv[3]
    bp = base_primes(int(math.isqrt(L1)) + 2)
    cnt = {}; last = None; ng = 0; t0 = time.time(); np_ = 0
    for s in range(L0, L1, 10 ** 8):
        P = odd_primes_in(s, min(s + 10 ** 8, L1), bp)
        np_ += len(P)
        if last is not None:
            P = np.concatenate([[last], P])
        last = int(P[-1]); g = np.diff(P); ng += len(g)
        v, c = np.unique(g[g >= 200], return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            cnt[a] = cnt.get(a, 0) + b
    print('primes', np_, 'gaps', ng, 'time', round(time.time() - t0), flush=True)
    json.dump({'L0': L0, 'L1': L1, 'primes': np_, 'gaps': ng, 'counts': {str(k): v for k, v in sorted(cnt.items())}},
              open(out, 'w'), indent=0)
