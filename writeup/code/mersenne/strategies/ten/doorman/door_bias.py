"""The doorman's slow mean, tested: a Chebyshev-type bias model.  For each non-trivial real character the sum over primes
<= x is about -1/2 pi(sqrt x) (the prime squares in the explicit formula; heuristic, needs x large against the modulus).
The set {p : n(p) > q_k} is described by k Legendre conditions, so its count is about M/2^k - 1/2 pi(sqrt x)(1 - 2^-k)
(M = number of odd primes up to x).  The model mean is 2 + sum_k (q_{k+1} - q_k) * tail_k / M with tail_k clipped at 0.
Compared with the measured tails and mean (ref21_*.json).   usage: door_bias.py"""
import json, math
from sympy import primepi, primerange
Q = list(primerange(2, 200))
E = sum(q / 2 ** (k + 1) for k, q in enumerate(Q))
for name, x in (('ref21_1e5.json', 10 ** 5), ('ref21_1e6.json', 10 ** 6)):
    R = json.load(open(name)); M = R['n']; cnt = {int(k): v for k, v in R['cnt'].items()}
    b = 0.5 * int(primepi(int(math.isqrt(x))))
    mean_meas = R['sum_n'] / M; mean_bias = 2.0
    print('x = %d, M = %d, 1/2 pi(sqrt x) = %.1f' % (x, M, b))
    print('   n(p) > q    measured tail   M/2^k    shortfall   bias model shortfall')
    for k in range(1, 12):
        q = Q[k - 1]; tail = sum(v for qq, v in cnt.items() if qq > q); exp = M / 2 ** k; model = max(0.0, exp - b * (1 - 2 ** -k))
        mean_bias += (Q[k] - q) * model / M
        if k <= 8: print('   %6d %14d %10.1f %10.1f %14.1f' % (q, tail, exp, exp - tail, b * (1 - 2 ** -k)))
    print('   mean: measured %.5f, bias model %.5f, Erdos %.5f; the model explains %.0f %% of the gap' % (mean_meas, mean_bias, E, 100 * (E - mean_bias) / (E - mean_meas)))
