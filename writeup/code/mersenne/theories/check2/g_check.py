# Item G, stage 3: check the pair-count machinery at large n by direct per-n counting,
# and print the four requested bins with more digits.
import numpy as np, os
d = os.path.dirname(os.path.abspath(__file__))
P = {a: np.load(os.path.join(d, 'G_p%d.npy' % a)) for a in (1, 5, 7, 11)}
bits = np.load(os.path.join(d, 'G_bits.npy'))
def isp_odd(x):
    i = x >> 1
    return (bits[i >> 3] >> (7 - (i & 7)).astype(np.uint8)) & 1
def G(a, b, t):
    pa = P[a][:np.searchsorted(P[a], t)]
    keys = (np.uint32(t - 1) - pa)[::-1]
    return int(np.searchsorted(P[b], keys, side='right').sum())
def direct(n, a):          # ordered pairs (p, n-p), p = a mod 12, both prime
    pa = P[a][P[a] < n].astype(np.int64)
    return int(isp_odd(n - pa).sum())
SPEC = {0: [(1, 1), (11, 1), (5, -1), (7, -1)], 6: [(1, 1), (5, 1), (7, -1), (11, -1)]}
lo = 501_187_234 + 1000
hi = lo + 12 * 15
for c, spec in SPEC.items():
    ns = [n for n in range(lo, hi) if n % 12 == c]
    dsum = sum(sg * direct(n, a) for n in ns for a, sg in spec)
    pairs = {0: [((1, 11), 2), ((5, 7), -2)], 6: [((1, 5), 2), ((7, 11), -2)]}[c]
    gsum = sum(w * (G(a, b, hi) - G(a, b, lo)) for (a, b), w in pairs)
    print('n = %d mod 12, %d n in [%d, %d): direct sum %d, via G %d, Delta values %s' % (
        c, len(ns), lo, hi, dsum, gsum, [sum(sg * direct(n, a) for a, sg in spec) for n in ns[:5]]))
r = np.load(os.path.join(d, 'G_results.npz')); e = r['edges']
for w in (1_000_000, 100_000_000, 446_683_592, 501_187_234, 891_250_938):
    i = int(np.nonzero(e == w)[0][0])
    print(w, ' '.join('%.4f' % r['D%d' % c][i] for c in (0, 2, 6, 10)), ' -S0 %.4f' % -r['S0'][i])
