# Line amplitudes A(rho) (dictionary normalisation: density term 2 pi A(rho) sqrt(rho) K(4 pi sqrt(rho y))),
# spike coefficients c(y0) of m(y) ~ c(y0)/sqrt(|y - y0|) (D<0: y > y0; D>0: y < y0), and a fine table of m_+-(y).
import os, sys, math, numpy as np
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import predlib as PL
z3 = PL.ZETA3
def s_cnt(v): return sum(1 for j in range(0, v // 3 + 1) if (v - 3 * j) % 2 == 0) if v >= 0 else 0
def beta(l, v):            # b_l(v)/b_l(0)
    if v == -1: return -1.0 / l**2
    if v == 0: return 1.0
    if v <= -2: return 0.0
    return (l - 1) / l * (s_cnt(v) + (1 - s_cnt(v - 1)) / l)
def fac(n):
    out = {}; q = 2
    while q * q <= n:
        while n % q == 0: out[q] = out.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out
def B(r):                  # B(rho) = A(rho) sqrt(rho)
    r = Fr(r); val = z3
    for l, v in fac(r.numerator).items(): val *= beta(l, v)
    for l, v in fac(r.denominator).items(): val *= beta(l, -v)
    return val
lines = []
for a in range(1, 2001):
    for b in range(1, 61):
        if math.gcd(a, b) == 1:
            Bv = B(Fr(a, b))
            if Bv != 0: lines.append((abs(Bv / math.sqrt(a / b)), Fr(a, b), Bv / math.sqrt(a / b), Bv))
lines.sort(key=lambda t: -t[0])
print('largest lines |A(rho)| (rho <= 2000, denominators <= 60):')
print('  ' + '; '.join('%s: %.4f' % (str(r), A) for _, r, A, _ in lines[:24]))
# spike coefficients: m ~ sum over (rho', k) of zeta(3) w(rho') rho'^(-1/2) / sqrt(|y - k^2/(4 rho')|)
def omega(a):
    val = 1.0
    for l, v in fac(a).items(): val *= PL.omega_l(l, v)
    return val
def wt(r):                 # w(rho')/zeta(3), rho' = a/b, b squarefree
    r = Fr(r); b = r.denominator; fb = fac(b)
    if any(v > 1 for v in fb.values()): return 0.0
    return omega(r.numerator) * np.prod([-1.0 / l**2 for l in fb]) if fb else omega(r.numerator)
spk = {}
alist = [(a, w) for (a, w, ps) in PL.enumerate_a(1e5, 1e-7)]
bl, mul = PL.squarefree_upto(200)
for a, w in alist:
    for b in bl:
        b = int(b)
        if math.gcd(a, b) != 1: continue
        wb = w * np.prod([-1.0 / l**2 for l in fac(b)]) if b > 1 else w
        rp = a / b; kmax = int(2 * math.sqrt(3 * rp))
        for k in range(1, kmax + 1):
            y0 = Fr(k * k * b, 4 * a)
            spk[y0] = spk.get(y0, 0.0) + z3 * wb / math.sqrt(rp)
big = sorted(spk.items(), key=lambda t: -abs(t[1]))[:16]
print('spike coefficients c(y0), y0 <= 3 (sums truncated: a <= 1e5 with |omega| a^-3/4 >= 1e-7, b <= 200):')
print('  ' + '; '.join('%s: %.4f' % (str(y0), c) for y0, c in sorted(big)))
# fine table: averages of m over bins of width 0.01
edges = np.round(np.linspace(0, 3, 301), 10)
with open(os.path.join(HERE, 'm_table.txt'), 'w') as f:
    f.write('# predicted murmuration density of S3 cubic fields, averages of m(y) = E[sqrt|D| a_K(p)] over bins of width 0.01 (limit X -> oo)\n')
    f.write('# y_lo y_hi m_minus(D<0) m_plus(D>0)\n')
    mm, _ = PL.predict_bins(edges, -1, Amax=3e5, tau=3e-6, B0=300, Bcap=5000)
    mp_, _ = PL.predict_bins(edges, 1, Amax=3e5, tau=3e-6, B0=300, Bcap=5000)
    for j in range(300): f.write('%.2f %.2f %9.4f %9.4f\n' % (edges[j], edges[j + 1], mm[j], mp_[j]))
print('m_table.txt written; samples (bin centre, m-, m+):', [(round((edges[j] + edges[j + 1]) / 2, 3), round(mm[j], 2), round(mp_[j], 2)) for j in (0, 4, 9, 24, 25, 49, 99, 100, 199, 224, 225, 299)])
