# Item B: Mertens predictor log-loss for lane words (c-1, c+1), c = 0 mod 6.
import numpy as np, math
from sympy import primerange, isprime

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return s

ISP = sieve(10_000_010)
small = list(primerange(5, 48))
mert = math.prod(1 - 1 / r for r in primerange(2, 48))
twin_tail = math.exp(sum(math.log1p(-1 / (r - 1) ** 2) for r in primerange(48, 2_000_000)))

def run(c_lo, c_hi, tail):
    c = np.arange(c_lo, c_hi + 1, 6, dtype=np.int64)
    assert np.all(c % 6 == 0)
    a, b = c - 1, c + 1
    cls = 2 * ISP[a].astype(int) + ISP[b].astype(int)
    sa = np.ones(len(c), bool); sb = np.ones(len(c), bool)
    for r in small:
        sa &= (a % r != 0); sb &= (b % r != 0)
    q = 1 / (np.log(c) * mert)
    pboth = q * q * tail
    P = np.zeros((len(c), 4))
    both = sa & sb; onlya = sa & ~sb; onlyb = ~sa & sb; none = ~sa & ~sb
    P[both, 3] = pboth[both]; P[both, 2] = (q - pboth)[both]; P[both, 1] = (q - pboth)[both]
    P[both, 0] = (1 - 2 * q + pboth)[both]
    P[onlya, 2] = q[onlya]; P[onlya, 0] = 1 - q[onlya]
    P[onlyb, 1] = q[onlyb]; P[onlyb, 0] = 1 - q[onlyb]
    P[none, 0] = 1
    pt = P[np.arange(len(c)), cls]
    assert np.all(pt > 0)
    # spot check classes with sympy
    for i in range(0, len(c), len(c) // 50):
        assert cls[i] == 2 * isprime(int(a[i])) + isprime(int(b[i]))
    return len(c), float(np.mean(-np.log2(pt)))

print('Mertens product r<=47: %.6f   twin tail prod_{47<r<2e6}: %.8f' % (mert, twin_tail))
print('c in [1000002, 1999998]: words %d, mean log-loss %.4f bits' % run(1_000_002, 1_999_998, twin_tail))
print('c in [9000000, 9999996]: words %d, mean log-loss %.4f bits' % run(9_000_000, 9_999_996, twin_tail))
twin_tail2 = math.exp(sum(math.log1p(-1 / (r - 1) ** 2) for r in primerange(48, 20_000_000)))
print('   (tail to 2e7: %.8f) -> %.4f' % (twin_tail2, run(9_000_000, 9_999_996, twin_tail2)[1]))
