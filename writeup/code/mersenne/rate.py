"""Hit rate after trial factoring: the chain p -> next prime q > p with no factor 2kq + 1 (k <= K) of 2^q - 1.
Every Mersenne exponent survives, so the chain visits all of them; hit rate = exponents / steps.
Prediction (Wagstaff-type heuristic): P(M_q prime | no factor below Q) ~ e^gamma log2(Q) / q, Q = 2Kq + 1."""
import sys, math, numpy as np
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243]
X, K0 = int(sys.argv[1]), int(sys.argv[2])
P, kmin = np.load(f'kmin_{X}_{K0}.npy')
E = set(MEXP)
eg = math.exp(0.5772156649)
print(f'primes 3 <= p < {X}: {len(P)}, exponents {sum(1 for m in MEXP if 2 < m < X)}')
print(f'{"K":>7} {"steps":>6} {"hits":>5} {"rate":>6} {"pred hits":>9}   first 10 / 20 / 30 steps of the chain from 1')
for K in [0, 1, 4, 10, 100, 1000, 10000, 100000]:
    if K > K0: break
    keep = (kmin == 0) | (kmin > K) if K else np.ones(len(P), bool)
    S = [2] + [int(p) for p in P[keep]]
    h = [s in E for s in S]
    pred = 1 + sum(min(1.0, eg * math.log2(2 * max(K, 1) * p + 1) / p) for p in P[keep])
    f = lambda n: f'{sum(h[:n])}/{n}'
    print(f'{K:>7} {len(S):>6} {sum(h):>5} {sum(h)/len(S):6.1%} {pred:9.1f}   {f(10)} {f(20)} {f(30)}   misses in the first 20: {[s for s in S[:20] if s not in E]}')
