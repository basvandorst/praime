"""Check Tracer.tr against PARI traces of Hecke matrices, and Tracer.sum_ap2 against trace(T_p^2)."""
import sys
sys.path.insert(0, sys.argv[1])
import cypari2
from hurwitz import hurwitz12
from trace import Tracer
from sympy import primerange
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
H = hurwitz12(4 * 10 ** 6)
bad = tot = 0
for N in (11, 23, 37, 43, 101):
    for k in (2, 4):
        pari('mf=mfinit([%d,%d],1)' % (N, k))
        T = Tracer(N, k, H)
        ns = [n for n in list(range(1, 61)) + [121, 169, 289, 361, 841, 961, 2 * 49, 3 * 25 * 7] if n % N]
        for n in ns:
            ref = int(pari('trace(mfheckemat(mf,%d))' % n))
            tot += 1
            if ref != T.tr(n):
                bad += 1; print("MISMATCH tr", N, k, n, ref, T.tr(n))
        T.basis()
        mx = 0
        for p in primerange(2, 200):
            if p == N: continue
            ref = int(pari('trace(mfheckemat(mf,%d)^2)' % p))
            v, err = T.sum_ap2(p); mx = max(mx, err); tot += 1
            if v != ref:
                bad += 1; print("MISMATCH ap2", N, k, p, ref, v)
        print("N=%d k=%d dim=%d S=%s max dev from integer %.2e" % (N, k, T.dim, T.S, mx), flush=True)
print("checks", tot, "mismatches", bad)
