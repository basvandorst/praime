"""Own sieve + own conditional MI estimator I(g_{n+1}; g_{n-1} | g_n) for primes in [L, L+W).
Controls: conditional shuffle (permute g_{n+1} within each value of g_n: keeps both pair tables exactly, true CMI 0),
and half samples. usage: cmi.py L W"""
import sys, math, json
import numpy as np

def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def primes_in(L, R, bp):
    # odd numbers in [L, R), L odd
    n = (R - L) // 2
    m = np.ones(n, bool)
    for p in bp[1:]:
        if p * p >= R: break
        st = max(p * p, ((L + p - 1) // p) * p)
        if st % 2 == 0: st += p
        m[(st - L) // 2::p] = False
    return L + 2 * np.nonzero(m)[0].astype(np.int64)

def H(c):
    c = np.asarray(c, float); c = c[c > 0]; n = c.sum(); p = c / n
    return -np.sum(p * np.log2(p)), (len(c) - 1) / (2 * n * math.log(2))

def cmi(a, b, c, B=200):
    k3 = np.bincount((a.astype(np.int64) * B + b) * B + c, minlength=B ** 3)
    k3 = k3.reshape(B, B, B)
    ab = k3.sum(2); bc = k3.sum(0); bb = ab.sum(0)
    out = {}
    hs = [H(ab.ravel()), H(bc.ravel()), H(bb), H(k3.ravel())]
    plug = hs[0][0] + hs[1][0] - hs[2][0] - hs[3][0]
    mm = plug + hs[0][1] + hs[1][1] - hs[2][1] - hs[3][1]
    return plug, mm, int(np.count_nonzero(k3))

if __name__ == '__main__':
    L, W = int(float(sys.argv[1])), int(float(sys.argv[2]))
    bp = base_primes(int(math.isqrt(L + W)) + 2)
    seg = 10 ** 7; gl = []; prev = None
    for s in range(L + 1 if L % 2 == 0 else L, L + W, seg):
        P = primes_in(s, min(s + seg, L + W), bp)
        if prev is not None: P = np.concatenate([[prev], P])
        gl.append((np.diff(P) // 2).astype(np.uint8)); prev = P[-1]
    g = np.concatenate(gl); del gl
    a, b, c = g[:-2], g[1:-1], g[2:]
    res = {'triples': len(a)}
    res['full'] = cmi(a, b, c)
    h = len(a) // 2
    res['half1'] = cmi(a[:h], b[:h], c[:h]); res['half2'] = cmi(a[h:], b[h:], c[h:])
    rng = np.random.default_rng(1)
    for r in range(2):
        cs = c.copy()
        for v in np.unique(b):
            idx = np.nonzero(b == v)[0]
            cs[idx] = rng.permutation(c[idx])
        res['shuffle%d' % r] = cmi(a, b, cs)
    print(json.dumps(res))
