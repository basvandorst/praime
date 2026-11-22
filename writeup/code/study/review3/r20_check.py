"""Round 20 check: pgf G(z) = sum_j T_j (z-1)^j / lam^j. ln G as a series in w=(z-1)/lam: c_k = factorial-cumulant sums.
Truncated at order K, evaluate at roots of unity, FFT -> P(k). Radius of convergence of ln G series = min |root| of Q(w)."""
import json, math
import numpy as np
E = json.load(open('/home/bas/scratch/study-review3/enum.json'))
r20 = {}
for f in ('r20_H30.json',):
    r20 = json.load(open('/home/bas/praime/writeup/data/study/' + f))
def logseries(T, n):
    # coefficients of ln(sum T_j w^j), T_0 = 1, up to w^n
    a = np.zeros(n + 1); a[:min(len(T), n + 1)] = T[:n + 1]
    c = np.zeros(n + 1)
    for k in range(1, n + 1):   # c' = Q'/Q  -> k c_k = k a_k - sum_{i=1}^{k-1} i c_i a_{k-i}
        c[k] = a[k] - sum(i * c[i] * a[k - i] for i in range(1, k)) / k
    return c
for H, lams in ((30, (16.6, 17.8, 23.0, 27.6, 34.5)), (60, (17.8, 23.0, 34.5))):
    T = np.array(E['H%d' % H]['T'])
    c = logseries(T, 12)
    print('H', H, 'cumulant sums c_k (k=1..8):', ['%.4g' % x for x in c[1:9]])
    roots = np.roots(T[::-1])
    rmin = min(abs(roots))
    for lam in lams:
        z = np.exp(2j * np.pi * np.arange(64) / 64)
        w = (z - 1) / lam
        exact = np.real(np.fft.fft(np.polyval(T[::-1], w))) / 64
        out = []
        for K in (3, 4, 5, 6, 8):
            G = np.exp(sum(c[k] * w ** k for k in range(1, K + 1)))
            P = np.real(np.fft.fft(G)) / 64
            out.append('K%d maxerr %.1e' % (K, np.max(np.abs(P[:12] - exact[:12]))))
        print('  lam %.1f  min|root w| = %.4f, 2/lam = %.4f -> ln-series %s on |z|=1;  G(-1) = %.4f; P0 exact %.5f;' % (
            lam, rmin, 2 / lam, 'DIVERGES' if 2 / lam > rmin else 'converges', np.polyval(T[::-1], -2 / lam), exact[0]), '; '.join(out))
