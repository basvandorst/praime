"""Part 4: ceilings on the touch rate.
(a) proved: a base-2 Fermat test never detects a composite M_p (2^(M_p - 1) = 1 mod M_p because p | 2^p - 2 = M_p - 1
    and 2 has order p mod M_p); checked on composite M_p below.
(b) heuristic: after a factor search that rules out every factor below Q, P(M_p prime) ~ e^g log2(Q)/p for every p
    (Bayes on the Wagstaff prior; the a factor cancels). So the posterior of a candidate reaches pi only when
    log2 Q >= pi p / e^g; table of the bits needed.
(c) calibration: sum of e^g log2(2Kp+1)/p over the TF survivors vs the true exponent count, p < 10^4 and p < 10^5."""
import math, numpy as np
eg = math.exp(0.5772156649015329)
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243]
for p in [11, 23, 29, 37, 41, 43, 47, 53]:
    M = (1 << p) - 1
    print(f'p={p}: M_p composite, 2^(M-1) mod M = {pow(2, M - 1, M)}, 3^(M-1) mod M = {pow(3, M - 1, M)}')
print('\nTF bits needed for posterior pi:   pi = 1%      10%      50%     (sqrt(M_p) = p/2 bits)')
for p in [10 ** 4, 10 ** 5, 10 ** 6, 136279841]:
    print(f'  p = {p:>10}: ' + '  '.join(f'{pi * p / eg:9.0f}' for pi in (0.01, 0.1, 0.5)) + f'   {p / 2:.0f}')
for bits in (64, 77, 80):
    p = 136279841
    print(f'  p = {p}, TF to 2^{bits}: posterior {eg * bits / p:.2e} = 1 in {p / (eg * bits):,.0f}')
for X, f, Ks in [(10000, 'kmin_10000_10000000.npy', [100, 10 ** 4, 10 ** 5, 10 ** 7]), (100000, 'kmin_100000_1000000.npy', [100, 10 ** 4, 10 ** 6])]:
    P, km = np.load(f)
    for K in Ks:
        s = [int(p) for p, k in zip(P, km) if k == 0 or k > K]
        pred = 1 + sum(min(1.0, eg * math.log2(2 * K * p + 1) / p) for p in s)
        mx = max(eg * math.log2(2 * K * p + 1) / p for p in s if p > 1000)
        print(f'X={X} K={K}: {len(s) + 1} candidates, predicted {pred:.1f} exponents, true {sum(1 for m in MEXP if m < X)}; '
              f'largest posterior among p > 1000: {mx:.2%}')
