"""Worked examples for theory 9b: single even n = 0 mod 12, their Goldbach ways by class, and the prediction."""
import numpy as np
N = 100_000_100
s = np.ones(N + 1, bool); s[:2] = False; s[4::2] = False
for i in range(3, int(N ** 0.5) + 1, 2):
    if s[i]: s[i * i::2 * i] = False
P = np.nonzero(s)[0]; P = P[P >= 5]
cls = {a: P[P % 12 == a] for a in (1, 5, 7, 11)}
sq = np.array([int(x) * int(x) for x in P if x * x <= N], np.int64)
def ways(n, a):
    p = cls[a][cls[a] < n - 4]
    return int(np.sum(s[n - p]))
def row(n):
    g111, g1111, g57, g75 = ways(n, 1), ways(n, 11), ways(n, 5), ways(n, 7)
    up, dn = g111 + g1111, g57 + g75
    S = int(np.sum(s[n - sq[sq < n]]))
    Q = int(np.sum(sq <= n))
    pi = {a: int(np.searchsorted(cls[a], n, side='right')) for a in cls}
    D = (pi[5] - pi[1]) + (pi[7] - pi[11])
    pred = -S * D / (Q / 2)
    return n, up, dn, up - dn, S, Q, D, pred
print('     n            ways via class 1   ways via 5 and 7   gap   S(n)   Q(n)   race lead D(n)   predicted gap')
for n in (1_000_008, 10_000_008, 100_000_008):
    n, up, dn, d, S, Q, D, pr = row(n)
    print(f'{n:>13,} {up:>16,} {dn:>18,} {d:>6} {S:>6} {Q:>6} {D:>16} {pr:>15.1f}')
print('\nneighbours n = 0 mod 12 just below 100,000,008:')
rs = [row(100_000_008 - 12 * k) for k in range(1000)]
for r in rs[:8]:
    print(f'{r[0]:>13,} {r[1]:>16,} {r[2]:>18,} {r[3]:>6} {r[4]:>6} {r[5]:>6} {r[6]:>16} {r[7]:>15.1f}')
for m in (10, 100, 1000):
    print(f'average of {m:>4} neighbours: gap {np.mean([r[3] for r in rs[:m]]):8.1f} (sd of that mean {np.std([r[3] for r in rs[:m]]) / np.sqrt(m):6.1f}), predicted {np.mean([r[7] for r in rs[:m]]):8.1f}, squares alone {-np.mean([r[4] for r in rs[:m]]):8.1f}')
