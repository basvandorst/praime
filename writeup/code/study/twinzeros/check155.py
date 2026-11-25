# Independent check of the selection rule of round 155 (plain enumeration, exact fractions for the real parts).
# For pairs (n, n+h): beta_h(chi)/S(h) = sum over the two members of the mean of conj(chi) over the admissible n mod q.
# Claim: = 2 mu(q) conj(chi(h)) / prod_{p|q}(p-2) if chi is even, q odd squarefree and coprime to h; else 0.
# Characters mod q are built from scratch: all homomorphisms (Z/q)^* -> roots of unity, by brute force on generators.
import cmath, math, itertools, sys
from sympy import factorint, primitive_root, totient
def chars(q):
    # all Dirichlet characters mod q as dicts n -> complex; via CRT on prime powers
    fac = factorint(q)
    comps = []
    for p, e in fac.items():
        pe = p ** e
        units = [a for a in range(pe) if math.gcd(a, pe) == 1]
        if p == 2:
            if e == 1: cs = [{1: 1}]
            elif e == 2: cs = [{1: 1, 3: 1}, {1: 1, 3: -1}]
            else:
                # (Z/2^e)^* = <-1> x <5>
                cs = []
                m = pe // 4
                for a in (0, 1):
                    for b in range(m):
                        d = {}
                        for s in (0, 1):
                            for t in range(m):
                                d[(pow(-1, s, pe) * pow(5, t, pe)) % pe] = (-1) ** (a * s) * cmath.exp(2j * cmath.pi * b * t / m)
                        cs.append(d)
        else:
            g = primitive_root(pe); phi = len(units); cs = []
            for b in range(phi):
                d = {}
                for t in range(phi): d[pow(g, t, pe)] = cmath.exp(2j * cmath.pi * b * t / phi)
                cs.append(d)
        comps.append((pe, cs))
    out = []
    for combo in itertools.product(*[c for _, c in comps]):
        d = {}
        for n in range(q):
            if math.gcd(n, q) == 1:
                v = 1
                for (pe, _), c in zip(comps, combo): v *= c[n % pe]
                d[n] = v
        out.append(d)
    return out
def primitive(d, q):
    for m in range(1, q):
        if q % m == 0:
            # induced from modulus m iff chi(a) = 1 for all a = 1 mod m
            if all(abs(d[a] - 1) < 1e-9 for a in d if (a - 1) % m == 0): return False
    return True
def mu(q):
    f = factorint(q)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)
nchk = nfail = 0; worst = 0.0; nprim = 0
for h in (2, 4, 6, 8, 10, 12, 30):
    for q in range(3, 61):
        adm = [n for n in range(q) if math.gcd(n, q) == 1 and math.gcd(n + h, q) == 1]
        for d in chars(q):
            if not primitive(d, q): continue
            if h == 2: nprim += 1
            if not adm: beta = 0
            else:
                beta = sum(d[n].conjugate() for n in adm) / len(adm) + sum(d[(n + h) % q].conjugate() for n in adm) / len(adm)
            even = abs(d[q - 1] - 1) < 1e-9
            if even and q % 2 == 1 and mu(q) != 0 and math.gcd(q, h) == 1:
                pred = 2 * mu(q) * d[h % q].conjugate() / math.prod(p - 2 for p in factorint(q))
            else: pred = 0
            err = abs(beta - pred); worst = max(worst, err)
            nchk += 1; nfail += (err > 1e-9)
print('primitive characters of conductor 3..60:', nprim, ' (h, chi) cases:', nchk, ' largest |enumeration - closed form|: %.2e' % worst)
# the inadmissible case: when no n is admissible (q = 2 for odd h) nothing is counted; q = 3, h = 2: one class
adm = [n for n in range(3) if n % 3 and (n + 2) % 3]
nchk += 1; nfail += (adm != [2]); print('twins mod 3: admissible classes', adm)
print('%d checks, %d failed' % (nchk, nfail))
