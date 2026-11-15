"""Exponents lit in the ord2 layout (prime order of 2 mod p) vs random exponent sets: 30 odd primes <= 500,000 drawn with the
heuristic weights w(p) = e^gamma log2(a p)/p (w = 1 for p <= 7), 4,000 sets; one-sided p-value."""
import numpy as np, math, sympy
N = 500000
E = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091]
pr = np.array(list(sympy.primerange(3, N + 1)))
lit = np.array([sympy.isprime(sympy.n_order(2, int(q))) for q in pr])
w = np.array([1.0 if p <= 7 else min(1.0, math.exp(0.5772156649) * math.log2((6 if p % 4 == 1 else 2) * p) / p) for p in pr])
idx = {int(p): i for i, p in enumerate(pr)}
real = int(sum(lit[idx[p]] for p in E))
rng = np.random.default_rng(11); sims = np.array([lit[rng.choice(len(pr), size=len(E), replace=False, p=w / w.sum())].sum() for _ in range(4000)])
print(f'lit exponents: {real} of {len(E)}; random sets: mean {sims.mean():.2f}, 95% {np.quantile(sims, .95):.0f}; P(>= {real}) = {(1 + (sims >= real).sum()) / 4001:.3f}')
