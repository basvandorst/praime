"""Question 3, method A: the classical criteria for non-idoneal P.

Representability: direct search over (a, b) (set of all a^2 + P b^2 <= N).
Criteria:
  P=27: p = 1 mod 3 and 2^((p-1)/3) = 1 mod p          (2 is a cubic residue)
  P=64: p = 1 mod 4 and 2^((p-1)/4) = 1 mod p          (2 is a 4th-power residue)
  P=14: (-14/p) = 1 and (x^2+1)^2 = 8 has a root: 8 = s^2 and s-1 or -s-1 a square
  P=11: x^3 - x^2 + x + 1 has three roots mod p        (cubic field of disc -44)
  P=23: x^3 - x - 1 has three roots mod p               (cubic field of disc -23)
Roots of the cubics: degree of gcd(f, x^p - x) with hand-written polynomial
arithmetic mod p.  Every prime p <= N is tested; primes dividing 4P or the
polynomial discriminant are listed separately.
Checked independently by qf_extra_check.py (PARI class polynomials, qfbsolve,
sympy galoistools).
"""
import json
from math import isqrt

import numpy as np

N = 500_000


def sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def represented(P):
    hit = np.zeros(N + 1, dtype=bool)
    b = 0
    while P * b * b <= N:
        a = np.arange(isqrt(N - P * b * b) + 1, dtype=np.int64)
        hit[a * a + P * b * b] = True
        b += 1
    return hit


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_mod_p(a, p):
    """Tonelli-Shanks; a a nonzero square mod odd prime p."""
    a %= p
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, bb * bb % p, t * bb * bb % p, r * bb % p
    return r


# --- polynomial arithmetic mod p (coefficient lists, lowest degree first) ---
def pmod(a, f, p):
    a = [x % p for x in a]
    df = len(f) - 1
    inv = pow(f[-1], p - 2, p)
    while len(a) - 1 >= df and any(a):
        if a[-1] == 0:
            a.pop()
            continue
        c = a[-1] * inv % p
        sh = len(a) - 1 - df
        for i in range(df + 1):
            a[sh + i] = (a[sh + i] - c * f[i]) % p
        a.pop()
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pmul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return r


def ppowx(e, f, p):
    """x^e mod (f, p)."""
    result, base = [1], [0, 1]
    while e:
        if e & 1:
            result = pmod(pmul(result, base, p), f, p)
        base = pmod(pmul(base, base, p), f, p)
        e >>= 1
    return result


def pgcd(a, b, p):
    a = pmod(a, [1], p) if False else [x % p for x in a]
    b = [x % p for x in b]
    while len(b) > 1 or (len(b) == 1 and b[0] != 0):
        while len(b) > 1 and b[-1] == 0:
            b.pop()
        if len(b) == 1 and b[0] == 0:
            break
        a, b = b, pmod(a, b, p)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def nroots(f, p):
    """Number of distinct roots of f in F_p (f monic, integer coefficients)."""
    fp = [c % p for c in f]
    xp = ppowx(p, fp, p)
    g = xp + [0] * max(0, 2 - len(xp))
    g[1] = (g[1] - 1) % p
    while len(g) > 1 and g[-1] == 0:
        g.pop()
    if len(g) == 1 and g[0] == 0:      # x^p - x = 0 mod f: f splits fully
        return len(f) - 1
    d = pgcd(fp, g, p)
    return len(d) - 1


def disc_cubic(f):
    c, b, a, _ = f      # x^3 + a x^2 + b x + c
    return a * a * b * b - 4 * b ** 3 - 4 * a ** 3 * c - 27 * c * c + 18 * a * b * c


def main():
    primes = [int(p) for p in sieve(N)]
    res = {}

    def tally(name, P, crit, bad_primes):
        rep = represented(P)
        agree = disagree = 0
        exc = []
        special = []
        for p in primes:
            if p in bad_primes:
                special.append((p, bool(rep[p]), crit(p)))
                continue
            if bool(rep[p]) == crit(p):
                agree += 1
            else:
                disagree += 1
                exc.append(p)
        nrep = int(sum(rep[p] for p in primes))
        res[name] = {"P": P, "checked": agree + disagree, "agree": agree,
                     "disagree": disagree, "exceptions": exc[:20],
                     "represented": nrep, "special": special}
        print(name, res[name], flush=True)

    tally("P=27 cubic residue", 27,
          lambda p: p % 3 == 1 and pow(2, (p - 1) // 3, p) == 1, set())
    tally("P=64 biquadratic residue", 64,
          lambda p: p % 4 == 1 and pow(2, (p - 1) // 4, p) == 1, set())

    def crit14(p):
        if legendre(-14, p) != 1:
            return False
        if legendre(8, p) != 1:
            return False
        s = sqrt_mod_p(8, p)
        return any(legendre(t, p) >= 0 for t in (s - 1, -s - 1))
    tally("P=14 (x^2+1)^2=8", 14, crit14, {2, 7})

    f11 = [1, 1, -1, 1]       # x^3 - x^2 + x + 1
    f23 = [-1, -1, 0, 1]      # x^3 - x - 1
    print("disc(x^3-x^2+x+1) =", disc_cubic(f11), " disc(x^3-x-1) =", disc_cubic(f23))
    tally("P=11 x^3-x^2+x+1 three roots", 11, lambda p: nroots(f11, p) == 3, {2, 11})
    tally("P=11 (-11/p)=1 and x^3-x^2+x+1 has a root", 11,
          lambda p: legendre(-11, p) == 1 and nroots(f11, p) >= 1, {2, 11})
    tally("P=23 x^3-x-1 three roots", 23, lambda p: nroots(f23, p) == 3, {2, 23})
    # the same criteria with polynomials instead of power residues (consistency)
    tally("P=27 x^3-2 three roots", 27, lambda p: nroots([-2, 0, 0, 1], p) == 3, {2, 3})
    tally("P=64 (p=1 mod 4) and x^4-2 has a root", 64,
          lambda p: p % 4 == 1 and nroots([-2, 0, 0, 0, 1], p) >= 1, {2})
    tally("P=14 x^4+2x^2-7 has a root and (-14/p)=1", 14,
          lambda p: legendre(-14, p) == 1 and nroots([-7, 0, 2, 0, 1], p) >= 1, {2, 7})
    # distribution of root counts of x^3-x-1 (Chebotarev in S3: 1/6, 1/2, 1/3)
    dist = {}
    for p in primes:
        if p != 23:
            k = nroots(f23, p)
            dist[k] = dist.get(k, 0) + 1
    res["x^3-x-1 root-count distribution (p != 23)"] = dist
    print("x^3-x-1 roots distribution", dist)
    json.dump(res, open("qf_extra.json", "w"), indent=1)


if __name__ == "__main__":
    main()
