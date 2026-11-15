"""Independent recheck helpers: no harness, no numpy sieve. Primes from sympy, weights from math,
exact Poisson-binomial by a plain-Python DP, expectations on the large ranges from stratified samples
of uniformly random primes (sympy.isprime) times exact prime counts (sympy.primepi)."""
import math, random
import sympy

EG = math.exp(0.5772156649015329)
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941,
        11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221,
        3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609,
        57885161, 74207281, 77232917, 82589933]
B = 83880077


def w(p):
    if p <= 7:
        return 1.0
    a = 6 if p % 4 == 1 else 2
    return min(1.0, EG * math.log2(a * p) / p)


def pb_tail(ws, O):
    """P(sum Bernoulli(ws) >= O), exact DP."""
    d = [1.0] + [0.0] * (O + 1)
    for x in ws:
        for j in range(O + 1, 0, -1):
            d[j] = d[j] * (1 - x) + d[j - 1] * x
        d[0] *= (1 - x)
    return 1.0 - sum(d[:O])


def poisson_tail(E, O):
    return 1.0 - sum(math.exp(-E) * E ** j / math.factorial(j) for j in range(O))


def stratified_E(lo, hi, member, per_bin=3000, seed=1):
    """sum of w(p)*member(p) over primes lo <= p < hi, estimated per doubling bin."""
    rnd = random.Random(seed)
    tot, var = 0.0, 0.0
    x = lo
    while x < hi:
        y = min(2 * x, hi)
        n = int(sympy.primepi(y - 1) - sympy.primepi(x - 1))
        vals = []
        while len(vals) < per_bin:
            c = rnd.randrange(x, y)
            if sympy.isprime(c):
                vals.append(w(c) * (1.0 if member(c) else 0.0))
        m = sum(vals) / len(vals)
        s2 = sum((v - m) ** 2 for v in vals) / (len(vals) - 1)
        tot += n * m
        var += n * n * s2 / len(vals)
        x = y
    return tot, math.sqrt(var)
