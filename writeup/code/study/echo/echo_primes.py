# The echo law, tested on primes alone (no zeros): with E(x) = psi(x) - x,
#   mean over ln x of E(a x) E(x) / x  =  a [1 + gamma + int_1^a E(t) dt/t^2] - ln 2pi + T(a)      (RH, simple zeros)
# Left side: primes up to XTOP, Hann-weighted mean over ln x in [ln X0, ln X1]. Right side: the prime powers up to a.
import sys, math, json, numpy as np
XTOP = float(sys.argv[1]); X0 = 1e3; AMAX = 12.0; X1 = XTOP / AMAX; DU = 1e-5; G = 0.5772156649015329
u = np.arange(math.log(X0), math.log(XTOP), DU); xs = np.exp(u); psi = np.zeros(len(u))
def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
B = base_primes(int(XTOP ** 0.5) + 2); SEG = 1 << 24; tot = 0.0; lo = 2
while lo < XTOP:
    hi = min(lo + SEG, int(XTOP) + 1); seg = np.ones(hi - lo, bool)
    for p in B.tolist():
        if p * p >= hi: break
        st = max(p * p, (lo + p - 1) // p * p); seg[st - lo::p] = False
    P = lo + np.nonzero(seg)[0]; c = np.concatenate(([0.0], np.cumsum(np.log(P))))
    i0, i1 = np.searchsorted(xs, lo), np.searchsorted(xs, hi)           # grid points with lo <= x < hi
    psi[i0:i1] = tot + c[np.searchsorted(P, xs[i0:i1], side='right')]; tot += c[-1]; lo = hi
pw = sorted((p ** k, math.log(p)) for p in B.tolist() for k in range(2, 40) if p ** k <= XTOP)   # prime powers
pn = np.array([t[0] for t in pw], float); pc = np.concatenate(([0.0], np.cumsum([t[1] for t in pw])))
psi += pc[np.searchsorted(pn, xs, side='right')]
E = psi - xs; n1 = int(round(math.log(X1 / X0) / DU)); w = np.hanning(n1); w /= w.sum(); f = E[:n1] / xs[:n1]
def T(a): return 1 - math.log(2) if a == 1 else 1 - 0.5 * math.log(1 - a ** -2) - a * math.atanh(1 / a)
def law(a):
    s = sum(math.log(p) * (1 / q - 1 / a) for p in (2, 3, 5, 7, 11) for q in [p ** k for k in range(1, 5)] if q <= a)
    return a * (1 + G + s - math.log(a)) - math.log(2 * math.pi) + T(a)
STEP = 100; ks = np.arange(0, int(math.log(AMAX) / DU), STEP); av = np.exp(ks * DU)
A = np.array([float(np.dot(w, E[k:k + n1] * f)) for k in ks]); L = np.array([law(a) for a in av])
C = 2 + G - math.log(4 * math.pi)
print('primes to %.1e, x from %.0e to %.1e (Hann mean over ln x, length %.1f); %d ratios a from 1 to %.0f' % (XTOP, X0, X1, math.log(X1 / X0), len(av), AMAX))
print('max |data - law| %.4f, rms %.4f; rms of the law itself %.4f' % (np.abs(A - L).max(), math.sqrt(np.mean((A - L) ** 2)), math.sqrt(np.mean(L ** 2))))
for a in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12):
    i = int(np.argmin(np.abs(av - a))); print('a = %5.2f  data %+.4f  law %+.4f  correlation %+.3f' % (av[i], A[i], L[i], L[i] / (math.sqrt(av[i]) * C)))
# corners: the jump of the slope at a prime power n is Lambda(n)/n; measured from straight-line fits 0.06 wide on each side
print('corner at n: jump of dA/da, data against Lambda(n)/n')
for n, lam in ((2, math.log(2)), (3, math.log(3)), (4, math.log(2)), (5, math.log(5)), (7, math.log(7)), (8, math.log(2)), (9, math.log(3)), (11, math.log(11)), (6, 0.0), (10, 0.0)):
    left = (av > n - 0.06) & (av < n - 0.002); right = (av > n + 0.002) & (av < n + 0.06)
    sl = np.polyfit(av[left], A[left], 1)[0]; sr = np.polyfit(av[right], A[right], 1)[0]
    print('  n = %2d  data %+.3f  law %+.3f' % (n, sr - sl, lam / n))
json.dump({'a': av.tolist(), 'data': A.tolist(), 'law': L.tolist()}, open('echo_%g.json' % XTOP, 'w'))
