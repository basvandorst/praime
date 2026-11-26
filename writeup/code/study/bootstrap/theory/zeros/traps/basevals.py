"""Base values for all degrees (numerics, 900-bit arb, own K0/K1):
m_j = g_j(0) (j <= 24), g_j(log n) for n = 2..NMAX, cost vector c_j = sum_n 2 Lambda(n) n^{-1/2} g_j(log n).
Output basevals.json (midpoints, 120 digits; 'acc_bits' = worst relative accuracy of the balls).
Usage: python basevals.py [NMAX=160]"""
import sys, json, time
from fam import *
from core import von_mangoldt
import mpmath as mp

J = 24
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 160
F = Fam(J, 0)
t0 = time.time()
r0 = F.rows(arb(1))
m = [r0[j][0] for j in range(J + 1)]
gl = {}
c = [arb(0)] * (J + 1)
last = None
for n in range(2, NMAX + 1):
    lam = von_mangoldt(n)
    r = F.rows(arb(n))
    gl[n] = [r[j][0] for j in range(J + 1)]
    if lam:
        w = 2 * arb(str(mp.nstr(lam, 5))) * 0   # placeholder, replaced below
        p = next(p for p in range(2, n + 1) if n % p == 0)
        w = 2 * arb(p).log() / arb(n).sqrt()
        for j in range(J + 1):
            c[j] = c[j] + w * gl[n][j]
        last = (n, max(abs(float((w * gl[n][j]).mid().str(5, radius=False)) if True else 0) for j in range(J + 1)))
print("last prime-power term (n, max_j |term_j|):", last)
print("max_j |c_j| =", max(float(abs(x).mid().str(5, radius=False)) for x in c))
acc = min(x.rel_accuracy_bits() for x in m + c)
print("min rel accuracy bits of m, c:", acc, " time %.0fs" % (time.time() - t0))
json.dump(dict(J=J, NMAX=NMAX, prec=ctx.prec, acc_bits=acc,
               m=[astr(x, 120) for x in m], c=[astr(x, 120) for x in c],
               gl={n: [astr(x, 120) for x in gl[n]] for n in range(2, 41)}),
          open("basevals.json", "w"))
