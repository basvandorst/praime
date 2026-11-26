"""Referee step 2: local mean at l | M of xi_Delta(l) * F_l(Delta) over p in Z_l^x, Delta = Q(s^2 Q - 4p).
Own l-adic code; xi computed from actual class numbers of a representative (Martin's definition), F_l from the
order class-number formula; enclosure for classes not resolved mod l^J. Then an end-to-end check over actual primes."""
import math
from fractions import Fraction as Fr
import cypari2
from sympy import primerange, factorint
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)

def H(n):
    return Fr(str(pari.qfbhclassno(n)))

def local_data(D, l):
    """D < 0 a discriminant. returns (e, chi): e = v_l(f0), chi = kronecker(D0, l)."""
    D0 = int(pari.coredisc(D)); f2 = D // D0; f = math.isqrt(f2); assert f * f == f2
    e = 0
    while f % l == 0:
        f //= l; e += 1
    return e, int(pari.kronecker(D0, l))

def F(e, chi, l):
    x = Fr(1, l)
    return sum(x ** j for j in range(e)) + x ** e / (1 - chi * x)

def xi_cls(D, l):
    if D % (l * l) == 0:
        return l * H(-D // (l * l)) / H(-D) - 1
    return int(pari.kronecker(D, l)) - 1

def e_new(s, l):
    return Fr(-l, l + 1) if s % l == 0 else Fr(-l * l, l * l - 1)

J = {2: 14, 3: 8, 5: 6}
res = []
for l in (2, 3, 5):
    mod = l ** J[l]
    for Qbase in (10103, 2 * 10103, 3 * 10103, 9973, 2 * 9973, 3 * 9973, 10111):
        if Qbase % l == 0:
            continue
        for s in range(0, 7):
            cache = {}; tot = Fr(0); und = 0; n = 0
            for b in range(1, mod):
                if b % l == 0:
                    continue
                n += 1
                d = (Qbase * (s * s * Qbase - 4 * b)) % mod
                v = 0; u = d
                if d == 0:
                    und += 1; continue
                while u % l == 0:
                    u //= l; v += 1
                if J[l] - v < (4 if l == 2 else 2):
                    und += 1; continue
                key = (v, u % (16 if l == 2 else l ** 3))
                if key not in cache:
                    # representative negative discriminant with the same class mod l^J (and = 0/1 mod 4)
                    rep = d - mod * 7919 * 4
                    while rep % 4 not in (0, 1) or rep >= 0:
                        rep -= mod
                    e, chi = local_data(rep, l)
                    cache[key] = xi_cls(rep, l) * F(e, chi, l)
                tot += cache[key]
            lo = tot / n - Fr(3 * und, n); hi = tot / n + Fr(3 * und, n)
            pred = e_new(s, l)
            res.append((l, Qbase, s, float(lo), float(hi), float(pred), lo <= pred <= hi))
bad = [r for r in res if not r[-1]]
print('local means: %d cases, %d outside enclosure; max width %.2e' % (len(res), len(bad), max(r[4] - r[3] for r in res)))
for r in res[:6] + [r for r in res if r[0] == 2][:7]:
    print('  l=%d Q=%d s=%d  [%.5f, %.5f]  pred %.5f' % r[:6])
for r in bad:
    print('  BAD', r)
