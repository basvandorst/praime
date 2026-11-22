"""I_k table with a numerically stable large-q part. H_q(m) = ln(q-1) + ((q-m)/(q-1)) phi(1/(q-m)), phi(e) = -(1-e)ln(1-e)."""
import math, sys
import numpy as np
from sympy import primerange
sys.path.insert(0, '.')
from r61_check import coupon_exact

def A(q, m):
    e = 1.0 / (q - m)
    return ((q - m) / (q - 1)) * (-(1 - e) * np.log1p(-e))

def big_part(qs, k):
    n = qs - 1
    r1 = np.exp(k * np.log1p(-1 / n)); r2 = np.exp(k * np.log1p(-2 / n))
    EU = (n - 1) * r1; VU = np.maximum((n - 1) * r1 + (n - 1) * (n - 2) * r2 - EU ** 2, 0)
    mu = n - EU
    h = np.maximum(0.05, np.minimum(np.sqrt(VU) + 1, (qs - mu) / 10)) * 0.5
    f0 = A(qs, mu); d2 = (A(qs, mu + h) - 2 * f0 + A(qs, mu - h)) / h ** 2
    return A(qs, np.ones_like(qs)) - (f0 + 0.5 * d2 * VU)

QX = 3000
small = [int(p) for p in primerange(3, QX)]
big = np.array(list(primerange(QX, 2 * 10 ** 6)), float)
res = {}
for k in (1, 2, 3, 5, 10, 30, 100, 1000, 10000):
    s = sum(coupon_exact(q, k) for q in small)
    tb = big_part(big, k)
    res[k] = s + tb.sum()
    print('k=%5d I_k %.8f nats %.5f bits | q in (2e5,2e6] part %.2e, last term*q^3/k %.3f' % (k, res[k], res[k] / math.log(2), tb[big > 2e5].sum(), tb[-1] * big[-1] ** 3 / k), flush=True)
print('I1 %.9f  I2-I1 %.9f nats' % (res[1], res[2] - res[1]))
print('growth 1e3->1e4 %.4f nats (lnln %.4f); 1e2->1e3 %.4f (lnln %.4f); 10->100 %.4f (lnln %.4f)' % (
    res[10000] - res[1000], math.log(math.log(1e4) / math.log(1e3)), res[1000] - res[100], math.log(math.log(1e3) / math.log(1e2)),
    res[100] - res[10], math.log(math.log(100) / math.log(10))))
