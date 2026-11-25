"""Own segmented odd-only sieve for the twin-class study (no repo code imported).

count_window(L0, W, hs, Q, nsub) counts, for n in [L0, L0 + W) split into nsub equal sub-windows:
  P[s, c]   number of primes p in sub-window s with p = c mod Q
  WP[s, c]  sum of 1/ln p over the same primes            (weight of the linear-response regressor)
  X[h][s, c] number of n = c mod Q in sub-window s with n and n + h both prime  (h even)
Primes p in [L0 + W, L0 + W + max h) are needed only as second members and are not put in P.
Usage as a script:  tcsieve.py L0 W out.npz [nsub]
"""
import sys, math, time
import numpy as np


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def count_window(L0, W, hs=(2, 4, 6, 8, 10, 12), Q=30030, nsub=9, seg_odd=1 << 23, verbose=True):
    L0 = int(L0); W = int(W)
    hmax = max(hs)
    look = hmax // 2                                   # extra odd entries needed after the segment
    hi = L0 + W + hmax
    R = small_primes(math.isqrt(hi) + 1)
    R = R[R > 2]
    sub_edges = [L0 + (W * s) // nsub for s in range(nsub + 1)]
    P = np.zeros((nsub, Q), np.int64); WP = np.zeros((nsub, Q)); X = {h: np.zeros((nsub, Q), np.int64) for h in hs}
    start = L0 + 1 if L0 % 2 == 0 else L0              # first odd >= L0
    t0 = time.time(); seg = start; nseg = 0
    Rl = R.tolist()
    while seg < L0 + W:
        n = min(seg_odd, (L0 + W - seg + 1) // 2)      # odd numbers seg, seg+2, ..., < L0+W
        ne = n + look
        a = np.ones(ne, dtype=bool)
        for r in Rl:
            if r * r > seg + 2 * ne:
                break
            j = ((-seg) * ((r + 1) // 2)) % r         # seg + 2j = 0 mod r
            if seg + 2 * j == r:                        # never strike r itself (only matters for tiny seg)
                j += r
            a[j::r] = False
        idx = np.nonzero(a[:n])[0]
        p = seg + 2 * idx
        s_of = np.searchsorted(np.array(sub_edges[1:]), p, side='right')
        cls = p % Q
        key = s_of * Q + cls
        P += np.bincount(key, minlength=nsub * Q).reshape(nsub, Q)
        WP += np.bincount(key, weights=1.0 / np.log(p.astype(np.float64)), minlength=nsub * Q).reshape(nsub, Q)
        for h in hs:
            m = a[idx + h // 2]
            X[h] += np.bincount(key[m], minlength=nsub * Q).reshape(nsub, Q)
        seg += 2 * n; nseg += 1
        if verbose and nseg % 50 == 0:
            print('  seg', nseg, 'at %.3e' % seg, '%.1fs' % (time.time() - t0), flush=True)
    if verbose:
        print('  done %.3e..%.3e in %.1fs: primes %d, twins %d' % (L0, L0 + W, time.time() - t0, P.sum(), X[hs[0]].sum()),
              flush=True)
    return P, WP, X


if __name__ == '__main__':
    L0, W, out = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    nsub = int(sys.argv[4]) if len(sys.argv) > 4 else 9
    P, WP, X = count_window(int(L0), int(W), nsub=nsub)
    np.savez_compressed(out, L0=int(L0), W=int(W), P=P.astype(np.int32), WP=WP,
                        **{'X%d' % h: X[h].astype(np.int32) for h in X})
