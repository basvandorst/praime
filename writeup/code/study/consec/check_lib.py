"""check_lib.py: independent checks of consec_lib labels (plain Python integers, sympy, PARI polclass / point counts).
Ranges: all odd primes 3 <= p < 3e4, every third prime in (1e8 - 1e5, 1e8) (dihedral, Gaussian), a_p for p < 6000."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sympy
from consec_lib import segment, pari

def roots_mod(coeffs, p):     # number of distinct roots of an integer polynomial mod p (PARI polrootsmod; independent of the forms)
    pol = 'Pol(%s)' % str(list(coeffs)).replace(' ', '')
    return int(pari()('#polrootsmod(%s, %d)' % (pol, p)))

def rep_brute(A, B, C, n):    # n = A x^2 + B x y + C y^2 with integers x, y (plain Python)
    D = B * B - 4 * A * C; ymax = math.isqrt(4 * A * n // (-D)) + 1
    for y in range(-ymax, ymax + 1):
        disc = B * B * y * y - 4 * A * (C * y * y - n)
        if disc < 0: continue
        s = math.isqrt(disc)
        if s * s != disc: continue
        for sg in (1, -1):
            if (-B * y + sg * s) % (2 * A) == 0: return True
    return False

P = pari(); bad = 0
H47 = P('polclass(-47)'); H39 = P('polclass(-39)')
for (L, U) in ((0, 30000), (10 ** 8 - 100000, 10 ** 8)):
    seg = segment(L, U, cases=('C23', 'C31', 'C47', 'C39', 'G'))
    p = seg['p']; ref = list(sympy.primerange(max(L, 3), U))
    assert list(p.tolist()) == ref, 'primes differ'
    sel = range(len(p)) if L == 0 else range(0, len(p), 3)
    for i in sel:
        q = int(p[i])
        # S3 cases: roots of the cubic
        for D, cub in ((23, [1, 0, -1, -1]), (31, [1, 0, 1, -1])):
            lab = int(seg['C%d' % D][i])
            if q == D: bad += lab != -1; continue
            want = {3: 1, 0: 2, 1: 0}[roots_mod(cub, q)]
            if lab != want: bad += 1; print('C%d' % D, q, lab, want)
        # D5: Hilbert class polynomial roots (principal) and brute-force forms
        if q != 47:
            lab = int(seg['C47'][i]); leg = sympy.legendre_symbol(q % 47, 47) if q % 47 else 0
            if leg == -1: want = 0
            else:
                nr = int(P('#polrootsmod(%s, %d)' % (H47, q)))
                want = 1 if nr == 5 else (2 if rep_brute(2, 1, 6, q) else 3)
                if want == 3: assert rep_brute(3, 1, 4, q) and nr == 0
            if lab != want: bad += 1; print('C47', q, lab, want)
        if q not in (3, 13):
            lab = int(seg['C39'][i]); chi = sympy.jacobi_symbol(q % 39, 39) if math.gcd(q, 39) == 1 else 0
            # (-39/q) = (q/3)(q/13) = jacobi(q, 39)
            if chi == -1: want = 0
            else:
                nr = int(P('#polrootsmod(%s, %d)' % (H39, q)))
                want = 1 if nr == 4 else (2 if rep_brute(3, 3, 4, q) else 3)
                if want == 3: assert rep_brute(2, 1, 5, q)
                if want == 1: assert rep_brute(1, 1, 10, q)
            if lab != want: bad += 1; print('C39', q, lab, want)
        # Gaussian
        th = seg['G_th'][i]
        if q % 4 == 3: bad += int(seg['G'][i]) != 0
        else:
            a = round(math.sqrt(q) * math.cos(th)); b = round(math.sqrt(q) * math.sin(th))
            ok = a * a + b * b == q and a > b > 0
            ok &= int(seg['G'][i]) == 1 + int(math.atan2(b, a) // (math.pi / 16))
            if not ok: bad += 1; print('G', q, th)
    print('segment', L, U, 'checked', len(sel), 'primes; bad so far', bad, flush=True)

# a_p of 37a1 against a point count
seg = segment(0, 6000, cases=('E37',))
for q, ap in zip(seg['p'].tolist(), seg['E37_ap'].tolist()):
    if q == 37: continue
    cnt = 1                                   # point at infinity
    for x in range(q):
        r = (x * x * x - x) % q; d = (1 + 4 * r) % q       # y^2 + y = r  <=>  (2y + 1)^2 = 1 + 4r
        cnt += 1 + (sympy.legendre_symbol(d, q) if d else 0)
    if q + 1 - cnt != ap: bad += 1; print('E37', q, ap, q + 1 - cnt)
print('a_p checked for', len(seg['p']), 'primes below 6000')
print('TOTAL BAD', bad)
