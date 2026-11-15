"""Question 3, method B (independent of qf_extra.py).

Representability: PARI qfbsolve(Qfb(1,0,P), p) (each solution re-verified).
Criterion 1 (theorem, Cox Thm 9.2 + Thm 11.1 / Cor. 11.37): for odd p not
  dividing P * disc(H), p = x^2 + P y^2 <=> (-P/p) = 1 and the class polynomial
  H_{-4P}(x) (PARI polclass, minimal polynomial of j(sqrt(-P))) has a root mod p.
  Roots counted with PARI polrootsmod.
Criterion 2: the simple polynomials of the classical statements, root counts by
  sympy.polys.galoistools (gcd with x^p - x), compared over all primes <= N.
Field check: PARI nfisisom(simple polynomial, polredabs(H)) nonempty, so both
  polynomials define the same field and the same splitting condition.
Primes from sympy.primerange.
"""
import json

import cypari2
from sympy import primerange, ZZ
from sympy.polys.galoistools import gf_pow_mod, gf_sub, gf_gcd, gf_from_int_poly

pari = cypari2.Pari()
pari.allocatemem(200 * 10**6)
N = 500_000
PRIMES = [int(p) for p in primerange(2, N + 1)]


def nroots_sympy(coeffs_high_first, p):
    f = gf_from_int_poly(coeffs_high_first, p)
    if len(f) - 1 < 1:
        return 0
    xp = gf_pow_mod([1, 0], p, f, p, ZZ)
    g = gf_sub(xp, [1, 0], p, ZZ)
    if not g:
        return len(f) - 1
    d = gf_gcd(f, g, p, ZZ)
    return len(d) - 1


def kron(a, p):
    return int(pari.kronecker(a, p))


SIMPLE = {
    27: ("p = 1 mod 3 and x^3 - 2 has a root", [1, 0, 0, -2],
         lambda p, r: p % 3 == 1 and r >= 1),
    64: ("p = 1 mod 4 and x^4 - 2 has a root", [1, 0, 0, 0, -2],
         lambda p, r: p % 4 == 1 and r >= 1),
    14: ("(-14/p) = 1 and (x^2+1)^2 - 8 has a root", [1, 0, 2, 0, -7],
         lambda p, r: kron(-14, p) == 1 and r >= 1),
    11: ("x^3 - x^2 + x + 1 has three roots", [1, -1, 1, 1],
         lambda p, r: r == 3),
    23: ("x^3 - x - 1 has three roots", [1, 0, -1, -1],
         lambda p, r: r == 3),
}
POLYSTR = {27: "x^3-2", 64: "x^4-2", 14: "x^4+2*x^2-7", 11: "x^3-x^2+x+1",
           23: "x^3-x-1"}


def main():
    out = {}
    for P in (27, 64, 14, 11, 23):
        Q = pari.Qfb(1, 0, P)
        rep = {}
        for p in PRIMES:
            s = pari.qfbsolve(Q, p)
            if len(s) > 0:
                x, y = int(s[0]), int(s[1])
                assert x * x + P * y * y == p
                rep[p] = True
            else:
                rep[p] = False
        H = pari.polclass(-4 * P)
        dH = pari.poldisc(H)
        Hbad = set(int(q) for q in pari.factor(abs(dH))[0]) | set(
            int(q) for q in pari.factor(4 * P)[0])
        a1 = d1 = 0
        for p in PRIMES:
            if p in Hbad:
                continue
            c = kron(-P, p) == 1 and len(pari.polrootsmod(H, p)) > 0
            if c == rep[p]:
                a1 += 1
            else:
                d1 += 1
        name, coeffs, crit = SIMPLE[P]
        a2 = d2 = 0
        exc2 = []
        for p in PRIMES:
            c = crit(p, nroots_sympy(coeffs, p))
            if c == rep[p]:
                a2 += 1
            else:
                d2 += 1
                exc2.append(p)
        iso = len(pari.nfisisom(POLYSTR[P], pari.polredabs(H))) > 0
        out[P] = {
            "represented": sum(rep.values()),
            "classpoly_degree": int(pari.poldegree(H)),
            "classpoly_excluded_primes": sorted(Hbad),
            "classpoly_checked": a1 + d1, "classpoly_agree": a1,
            "classpoly_disagree": d1,
            "simple_criterion": name, "simple_checked": a2 + d2,
            "simple_agree": a2, "simple_disagree": d2,
            "simple_exceptions": exc2,
            "polredabs(H)": str(pari.polredabs(H)),
            "same_field_as_simple": iso,
        }
        print(P, out[P], flush=True)
    json.dump(out, open("qf_extra_check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
